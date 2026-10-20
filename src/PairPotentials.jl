module PairPotentials
using SIMD: vifelse, Vec
using AtomsBase: ChemicalSpecies
using ..NeighbourList: PairBlock, pack
export AbstractPairPotential, EnergyOnly, EnergyForcesVirial, SmoothPair, RadialPotential, cutoff, bind_species

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
@inline select(mask::Bool,a,b) = ifelse(mask,a,b)
@inline select(mask::Vec,a,b) = vifelse(mask,a,b)

"""Smooth finite-cutoff example: U_ab(r)=A_ab*(1-(r/c)^2)^3 for r<c.
`labels` are AtomsBase species values; symmetric amplitudes in eV, cutoff in Å.
No implicit mixing rule. Values and species are runtime data.
"""
struct SmoothPair{T} <: AbstractPairPotential
    labels::Vector{ChemicalSpecies}
    amplitudes::Matrix{T}
    radius::T
    SmoothPair{T}(labels, amplitudes, radius) where T = new{T}(labels, amplitudes, radius)
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
struct BoundSmooth{T}
    amplitudes::Matrix{T}
    radius::T
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
end
