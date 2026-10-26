# Defined in parameters.jl, included after this calculator boundary.
function parameter_energy end

AC.energy_unit(::AbstractPairPotential)=u"eV"
AC.length_unit(::AbstractPairPotential)=u"Å"

# Select precision only at the system boundary; geometry conversion stays explicit
# in prepare!. Empty systems do not require position(sys,1).
function default_precision(sys)
    isempty_system = length(sys)==0
    isempty_system && return Float64
    typeof(ustrip(position(sys,1)[1]))===Float32 ? Float32 : Float64
end
function AC.promote_force_type(sys::AbstractSystem{3}, ::AbstractPairPotential)
    SVector{3,typeof(zero(default_precision(sys))*u"eV/Å")}
end
AC.zero_energy(sys::AbstractSystem{3}, ::AbstractPairPotential)=zero(default_precision(sys))*u"eV"
function AC.zero_virial(sys::AbstractSystem{3}, p::AbstractPairPotential)
    zero(SMatrix{3,3,typeof(AC.zero_energy(sys,p)),9})
end

workspace_width(::Workspace{W}) where W = W
function select_workspace(sys, neighbourlist, width, precision, bruteforce_threshold=nothing)
    if isnothing(neighbourlist)
        return Workspace(width=something(width,4),precision=something(precision,default_precision(sys)),
            bruteforce_threshold=something(bruteforce_threshold,32))
    end
    neighbourlist isa Workspace || throw(ArgumentError("neighbourlist must be a NeighbourList.Workspace or nothing"))
    isnothing(width) || width==workspace_width(neighbourlist) ||
        throw(ArgumentError("width conflicts with the supplied neighbourlist workspace"))
    isnothing(precision) || precision===eltype(neighbourlist.x) ||
        throw(ArgumentError("precision conflicts with the supplied neighbourlist workspace"))
    isnothing(bruteforce_threshold) || bruteforce_threshold==neighbourlist.bruteforce_threshold ||
        throw(ArgumentError("bruteforce_threshold conflicts with the supplied neighbourlist workspace"))
    neighbourlist
end

function compute(sys,p,request; neighbourlist=nothing,width=nothing,precision=nothing,
                 parallel=false,parameters=nothing,bruteforce_threshold=nothing,kwargs...)
    (isnothing(parameters) || parameters===NamedTuple()) ||
        throw(ArgumentError("explicit parameters require calculate(Energy(), ...; parameters)"))
    # Unknown keywords are accepted for AtomsCalculators composition/conformance.
    w=select_workspace(sys,neighbourlist,width,precision,bruteforce_threshold)
    evaluate!(w,sys,p;request,parallel)
end

# Vector-of-static-vectors is the AtomsCalculators boundary representation.
# This conversion follows the existing one-time original-order permutation.
function force_vectors(f::AbstractMatrix{T}) where T
    [SVector{3,T}(f[1,i],f[2,i],f[3,i]) for i in axes(f,2)]
end
static_virial(v)=SMatrix{3,3}(v)

function AC.potential_energy(sys::AbstractSystem{3},p::AbstractPairPotential;kwargs...)
    compute(sys,p,EnergyOnly();kwargs...).energy
end
function AC.forces(sys::AbstractSystem{3},p::AbstractPairPotential;kwargs...)
    force_vectors(compute(sys,p,EnergyForcesVirial();kwargs...).forces)
end
function AC.virial(sys::AbstractSystem{3},p::AbstractPairPotential;kwargs...)
    static_virial(compute(sys,p,EnergyForcesVirial();kwargs...).virial)
end
function AC.energy_forces(sys::AbstractSystem{3},p::AbstractPairPotential;kwargs...)
    r=compute(sys,p,EnergyForcesVirial();kwargs...)
    (;energy=r.energy,forces=force_vectors(r.forces))
end
function AC.energy_forces_virial(sys::AbstractSystem{3},p::AbstractPairPotential;kwargs...)
    r=compute(sys,p,EnergyForcesVirial();kwargs...)
    (;energy=r.energy,forces=force_vectors(r.forces),virial=static_virial(r.virial))
end

function check_destination(f,sys)
    axes(f)==(Base.OneTo(length(sys)),) || throw(DimensionMismatch("force destination must have one entry per atom"))
end
function add_forces!(f,increment)
    for i in eachindex(f,increment); f[i]+=increment[i]; end
    f
end
function AC.forces!(f::AbstractVector,sys::AbstractSystem{3},p::AbstractPairPotential;kwargs...)
    check_destination(f,sys)
    add_forces!(f,AC.forces(sys,p;kwargs...))
end
function AC.energy_forces!(f::AbstractVector,sys::AbstractSystem{3},p::AbstractPairPotential;kwargs...)
    check_destination(f,sys)
    r=AC.energy_forces(sys,p;kwargs...)
    add_forces!(f,r.forces)
    (;energy=r.energy,forces=f)
end
function AC.energy_forces_virial!(f::AbstractVector,sys::AbstractSystem{3},p::AbstractPairPotential;kwargs...)
    check_destination(f,sys)
    r=AC.energy_forces_virial(sys,p;kwargs...)
    add_forces!(f,r.forces)
    (;energy=r.energy,forces=f,virial=r.virial)
end

property_result(::AC.Energy,r)=(;energy=r.energy)
property_result(::AC.Forces,r)=(;forces=force_vectors(r.forces))
property_result(::AC.Virial,r)=(;virial=static_virial(r.virial))
const SingleProperty=Union{AC.Energy,AC.Forces,AC.Virial}
function calculate(properties,sys,p,parameters,state;kwargs...)
    props=properties isa Tuple ? properties : (properties,)
    normalized=NL.normalize_properties(props)
    all(prop -> prop isa Union{AC.Energy,AC.Forces,AC.Virial},props) ||
        throw(ArgumentError("supported properties are Energy(), Forces() and Virial()"))
    if !(isnothing(parameters) || parameters===NamedTuple())
        parameters isa NamedTuple || throw(ArgumentError("explicit parameters must be a NamedTuple"))
        isempty(props) && return (;state)
        if all(prop -> prop isa AC.Energy,props)
            energy=parameter_energy(sys,p,parameters;kwargs...)*u"eV"
            return (;energy,state)
        end
        # Explicit scalar correctness route. Production bundled-parameter SIMD
        # execution is unchanged. Never convert explicit parameters to workspace
        # precision or detach force/virial parameter arithmetic.
        input=pair_inputs(sys,p)
        r=pair_evaluate(input.positions,p,parameters,input.environment,input.cell)
        unitful=(;energy=r.energy*u"eV",forces=r.forces*u"eV/Å",virial=r.virial*u"eV")
        return merge(map(prop -> property_result(prop,unitful),props)...,(;state))
    end
    isempty(props) && return (;state)
    request=normalized
    r=compute(sys,p,request;kwargs...)
    merge(map(prop -> property_result(prop,r),props)...,(;state))
end

# Support both the documented positional parameters/state form and the keyword
# form used by AtomsCalculators.Testing. Tuple requests use one fused traversal.
for PropertyType in (SingleProperty,Tuple)
    @eval begin
        function AC.calculate(properties::$PropertyType,sys::AbstractSystem{3},p::AbstractPairPotential;
                              parameters=nothing,state=nothing,kwargs...)
            calculate(properties,sys,p,parameters,state;kwargs...)
        end
        function AC.calculate(properties::$PropertyType,sys::AbstractSystem{3},p::AbstractPairPotential,
                              parameters;state=nothing,kwargs...)
            calculate(properties,sys,p,parameters,state;kwargs...)
        end
        function AC.calculate(properties::$PropertyType,sys::AbstractSystem{3},p::AbstractPairPotential,
                              parameters,state;kwargs...)
            calculate(properties,sys,p,parameters,state;kwargs...)
        end
    end
end
