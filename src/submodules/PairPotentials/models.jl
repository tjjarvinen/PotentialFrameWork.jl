"""`with_parameters(potential, parameters::NamedTuple)` returns a new callback.
Extend together with `AtomsCalculators.get_parameters` for parameterised models.
Parameters are numeric (eV/Å conventions); the physical cutoff must stay fixed.
The reconstructed callback must accept scalar parameter-promoted PairBlocks.
"""
function with_parameters(p, parameters)
    throw(ArgumentError("$(typeof(p)) does not implement with_parameters for $(typeof(parameters))"))
end

"""Subtype to give a concrete pair-block callback the AtomsCalculators interface.
Implement `cutoff` and `(potential)(block, request)`; optionally `bind_species`.
Parameters are bundled in the potential; callbacks use numeric Å and eV.
"""
abstract type AbstractPairPotential end

struct EnergyOnly end
struct EnergyForcesVirial end
"""Physical cutoff in Å. Extend for each concrete interaction."""
function cutoff end
"""Bind runtime system species to potential parameters once, outside hot loops.
Custom interactions receive ChemicalSpecies values and may use the default."""
bind_species(p, species_map) = p

"""Smooth finite-cutoff example: U_ab(r)=A_ab*(1-(r/c)^2)^3 for r<c.
`labels` are AtomsBase species values; symmetric amplitudes in eV, cutoff in Å.
No implicit mixing rule. Values and species are runtime data.
"""
struct SmoothPair{T,R} <: AbstractPairPotential
    labels::Vector{ChemicalSpecies}
    amplitudes::Matrix{T}
    radius::R
    SmoothPair{T}(labels, amplitudes, radius::R) where {T,R} = new{T,R}(labels, amplitudes, radius)
end
function SmoothPair(labels, amplitudes::AbstractMatrix{T}, radius::Real) where {T<:AbstractFloat}
    labels=ChemicalSpecies.(labels)
    n=length(labels)
    size(amplitudes)==(n,n) || throw(ArgumentError("parameter table size must match species labels"))
    length(unique(labels))==n || throw(ArgumentError("species labels must be unique"))
    all(isfinite,amplitudes) && amplitudes==transpose(amplitudes) ||
        throw(ArgumentError("amplitude table must be finite and symmetric"))
    isfinite(radius) && radius>0 || throw(ArgumentError("cutoff must be finite and positive"))
    SmoothPair{T}(collect(labels),Matrix(amplitudes),T(radius))
end
cutoff(p::SmoothPair)=p.radius
struct BoundSmooth{T,R}
    amplitudes::Matrix{T}
    radius::R
    indices::Dict{ChemicalSpecies,Int}
end
function bind_species(p::SmoothPair, labels)
    ids=map(labels) do s
        id=findfirst(isequal(s),p.labels)
        isnothing(id) && throw(ArgumentError("potential has no parameters for species $s"))
        id
    end
    BoundSmooth(p.amplitudes[ids,ids],p.radius,
        Dict{ChemicalSpecies,Int}(s=>i for (i,s) in enumerate(labels)))
end
@inline function (p::BoundSmooth)(b::PairBlock{W,V}, req) where {W,V}
    # Convert parameters to geometry precision, keeping Float32 arithmetic Float32.
    T=eltype(V)
    central=p.indices[b.central_species]
    a=pack(Val(W), k -> T(p.amplitudes[central,p.indices[b.neighbour_species[k]]]))
    invc2=inv(T(p.radius)^2)
    q=1-(b.dx*b.dx+b.dy*b.dy+b.dz*b.dz)*invc2
    e=a*q*q*q
    req isa EnergyOnly && return (;energy=e)
    g=-6*a*q*q*invc2
    (;energy=e, fx=g*b.dx, fy=g*b.dy, fz=g*b.dz)
end

"""Adapter for `radial(r, central_species, neighbour_species, request)`.
For EnergyOnly return `(energy=U,)`; otherwise return `(energy=U, derivative=dUdr)`.
`r` is scalar for W=1 and SIMD.Vec for W>1. Units are Å/eV; masking is executor-owned.
At r=0 the adapter selects a safe denominator; radial models must have a regular
origin with dU/dr=0 there (or reject coincident interacting atoms).
"""
struct RadialPotential{R,T} <: AbstractPairPotential
    radial::R
    radius::T
end
cutoff(p::RadialPotential)=p.radius
@inline function (p::RadialPotential)(b::PairBlock, req)
    r=sqrt(b.dx*b.dx+b.dy*b.dy+b.dz*b.dz)
    out=p.radial(r,b.central_species,b.neighbour_species,req)
    req isa EnergyOnly && return (;energy=out.energy)
    safe=select(r>0,r,one(r))
    g=out.derivative/safe
    (;energy=out.energy, fx=g*b.dx, fy=g*b.dy, fz=g*b.dz)
end

# Legacy fx/fy/fz already represent dE/dd (the center force contribution).
# Energy-only avoids derivatives; other subsets may calculate/discard extras.
struct LegacyPair{P}
    potential::P
end
AC.length_unit(p::LegacyPair)=applicable(AC.length_unit,p.potential) ? AC.length_unit(p.potential) : u"Å"
AC.energy_unit(p::LegacyPair)=applicable(AC.energy_unit,p.potential) ? AC.energy_unit(p.potential) : u"eV"
NI.evaluation_kind(::AbstractPairPotential)=:pair
NI.search_cutoff(p::AbstractPairPotential)=cutoff(p)
NI.species_channels(::AbstractPairPotential)=nothing
function NI.prepare_potential(p::AbstractPairPotential,species;precision=Float64,width=1,kwargs...)
    NL.check_units(p)
    LegacyPair(bind_species(p,species))
end
NI.make_scratch(::LegacyPair,capacity,precision)=nothing
@inline function NI.evaluate_pair!(out,p::LegacyPair,block,ps,scratch,props)
    bound=isnothing(ps) ? p.potential : with_parameters(p.potential,ps)
    derivatives=NL.wants(AC.Forces,props) || NL.wants(AC.Virial,props)
    value=bound(block,derivatives ? EnergyForcesVirial() : EnergyOnly())
    NL.wants(AC.Energy,props) && NI.set_energy!(out,value.energy)
    derivatives && NI.set_gradient!(out,(value.fx,value.fy,value.fz))
    nothing
end

# Native radial adapter: expose only dE/dr to the executor. Keep the historical
# callable Cartesian adapter above for the functional AD and legacy entrypoints.
struct PreparedRadial{P}
    potential::P
end
AC.length_unit(p::PreparedRadial)=AC.length_unit(p.potential)
AC.energy_unit(p::PreparedRadial)=AC.energy_unit(p.potential)
function NI.prepare_potential(p::RadialPotential,species;kwargs...)
    NL.check_units(p)
    PreparedRadial(p)
end
NI.make_scratch(::PreparedRadial,capacity,precision)=nothing
@inline function NI.evaluate_pair!(out,p::PreparedRadial,block,ps,scratch,props)
    bound=isnothing(ps) ? p.potential : with_parameters(p.potential,ps)
    derivatives=NL.wants(AC.Forces,props) || NL.wants(AC.Virial,props)
    value=bound.radial(NI.distance(block),NI.central_species(block),NI.neighbour_species(block),
        derivatives ? EnergyForcesVirial() : EnergyOnly())
    NL.wants(AC.Energy,props) && NI.set_energy!(out,value.energy)
    derivatives && NI.set_radial_derivative!(out,value.derivative)
    nothing
end

struct PreparedSmooth{P}
    potential::P
    parameter_indices::Matrix{Int}
    parameter_count::Int
end
AC.length_unit(::PreparedSmooth)=u"Å"
AC.energy_unit(::PreparedSmooth)=u"eV"
function NI.prepare_potential(p::SmoothPair,species;kwargs...)
    NL.check_units(p)
    bound=bind_species(p,species)
    ids=[findfirst(isequal(s),p.labels) for s in species]
    parameter_indices=[let (a,b)=minmax(ids[i],ids[j]); b*(b-1)÷2+a; end for i in eachindex(ids),j in eachindex(ids)]
    PreparedSmooth(bound,parameter_indices,length(p.labels)*(length(p.labels)+1)÷2)
end
NI.make_scratch(::PreparedSmooth,capacity,precision)=nothing
@inline function NI.evaluate_pair!(out,prepared::PreparedSmooth,block,ps,scratch,props)
    p=prepared.potential; r=NI.distance(block); T=eltype(typeof(r))
    central=p.indices[NI.central_species(block)]
    species=NI.neighbour_species(block)
    a=pack(Val(NI.neighbour_count(block)),k -> begin
        neighbour=p.indices[species[k]]
        T(isnothing(ps) ? p.amplitudes[central,neighbour] : ps.amplitudes[prepared.parameter_indices[central,neighbour]])
    end)
    invc2=inv(T(p.radius)^2); q=1-r*r*invc2
    NL.wants(AC.Energy,props) && NI.set_energy!(out,a*q*q*q)
    if NL.wants(AC.Forces,props) || NL.wants(AC.Virial,props)
        NI.set_radial_derivative!(out,-6*a*q*q*invc2*r)
    end
    nothing
end

function NI.update_prepared!(p::PreparedSmooth,ps)
    isnothing(ps) && return nothing
    keys(ps)==(:amplitudes,) || throw(ArgumentError("SmoothPair parameters must contain only amplitudes; cutoff is fixed"))
    ps.amplitudes isa AbstractVector{<:Real} && length(ps.amplitudes)==p.parameter_count ||
        throw(ArgumentError("incorrect SmoothPair amplitude count"))
    all(isfinite,ps.amplitudes) || throw(ArgumentError("amplitudes must be finite"))
    nothing
end
