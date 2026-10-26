"""Convenient owning inference wrapper. `parameters` is authoritative; update it
with `AtomsCalculators.set_parameters!`. Environment is rebuilt per public call.
Do not share a calculator's mutable environment between concurrent calls.
"""
mutable struct ACECalculator{M,P,S,W}
    model::M
    parameters::P
    state::S
    environment::FullEnvironment
    executor::W
end
function ACECalculator(model,ps,st=NamedTuple();width=4,slots=Threads.nthreads(),precision=Float64,angular_backend=:sphericart,batch_size=1)
    env=FullEnvironment(;width)
    executor=CPUExecutor(model,env;slots,precision,angular_backend,batch_size)
    ACECalculator(model,deepcopy(ps),deepcopy(st),env,executor)
end
AC.get_parameters(c::ACECalculator)=deepcopy(c.parameters)
function AC.set_parameters!(c::ACECalculator,ps)
    size(ps.WB)==size(c.parameters.WB) && size(ps.Wpair)==size(c.parameters.Wpair) || throw(DimensionMismatch("ACE readout shapes"))
    c.parameters=deepcopy(ps)
    c
end
AC.energy_unit(::ACECalculator)=u"eV"
AC.length_unit(::ACECalculator)=u"Å"
AC.promote_force_type(sys::AbstractSystem{3},c::ACECalculator)=SVector{3,typeof(zero(eltype(c.executor.force))*u"eV/Å")}
AC.zero_energy(sys::AbstractSystem{3},c::ACECalculator)=zero(eltype(c.executor.force))*u"eV"
AC.zero_virial(sys::AbstractSystem{3},c::ACECalculator)=zero(SMatrix{3,3,typeof(AC.zero_energy(sys,c)),9})
function calculator_result(sys,c;parameters=c.parameters,energy_only=false,parallel=false,simd=true,fused=false,centers=nothing,properties=(AC.Energy(),AC.Forces(),AC.Virial()),kwargs...)
    parameters=(isnothing(parameters) || parameters===NamedTuple()) ? c.parameters : parameters
    fast_parameters=parameter_seed(parameters) isa Union{Float32,Float64}
    input=ace_inputs(sys,c.model;environment=c.environment)
    y=if energy_only || !fast_parameters
        ace_evaluate(input.positions,c.model,parameters,input.environment,input.cell;state=c.state,energy_only,centers)
    else
        prepare_executor!(c.executor,c.environment)
        execute_prepared!(c.executor,input.positions,input.cell,parameters;parallel,simd,fused,centers,properties)
        prepared_result(c.executor)
    end
    if energy_only
        e=fast_parameters ? eltype(c.executor.force)(y.energy) : y.energy
        return (;energy=e*u"eV")
    end
    (;energy=y.energy*u"eV",forces=force_vectors(y.forces*u"eV/Å"),virial=static_virial(y.virial*u"eV"))
end
AC.potential_energy(sys::AbstractSystem{3},c::ACECalculator;kwargs...)=calculator_result(sys,c;energy_only=true,kwargs...).energy
AC.forces(sys::AbstractSystem{3},c::ACECalculator;kwargs...)=calculator_result(sys,c;kwargs...).forces
AC.virial(sys::AbstractSystem{3},c::ACECalculator;kwargs...)=calculator_result(sys,c;kwargs...).virial
AC.energy_forces_virial(sys::AbstractSystem{3},c::ACECalculator;kwargs...)=calculator_result(sys,c;kwargs...)
function AC.energy_forces(sys::AbstractSystem{3},c::ACECalculator;kwargs...)
    r=calculator_result(sys,c;kwargs...); (;energy=r.energy,forces=r.forces)
end
function AC.forces!(f::AbstractVector,sys::AbstractSystem{3},c::ACECalculator;kwargs...)
    check_destination(f,sys); add_forces!(f,AC.forces(sys,c;kwargs...))
end
function AC.energy_forces!(f::AbstractVector,sys::AbstractSystem{3},c::ACECalculator;kwargs...)
    check_destination(f,sys); r=AC.energy_forces(sys,c;kwargs...)
    add_forces!(f,r.forces); (;energy=r.energy,forces=f)
end
function AC.energy_forces_virial!(f::AbstractVector,sys::AbstractSystem{3},c::ACECalculator;kwargs...)
    check_destination(f,sys); r=AC.energy_forces_virial(sys,c;kwargs...)
    add_forces!(f,r.forces); (;energy=r.energy,forces=f,virial=r.virial)
end
function calculate_ace(properties,sys::AbstractSystem{3},c::ACECalculator;
                      parameters=c.parameters,state=c.state,kwargs...)
    props=properties isa Tuple ? properties : (properties,)
    normalized=NL.normalize_properties(props)
    all(p -> p isa Union{AC.Energy,AC.Forces,AC.Virial},props) || throw(ArgumentError("unsupported ACE property"))
    isempty(props) && return (;state)
    r=calculator_result(sys,c;parameters,energy_only=all(p -> p isa AC.Energy,props),properties=normalized,kwargs...)
    key(::AC.Energy)=:energy; key(::AC.Forces)=:forces; key(::AC.Virial)=:virial
    merge((NamedTuple{(key(p),)}((r[key(p)],)) for p in props)...,(;state))
end
for PropertyType in (Union{AC.Energy,AC.Forces,AC.Virial},Tuple)
    @eval begin
        AC.calculate(props::$PropertyType,sys::AbstractSystem{3},c::ACECalculator;kwargs...)=calculate_ace(props,sys,c;kwargs...)
        AC.calculate(props::$PropertyType,sys::AbstractSystem{3},c::ACECalculator,parameters;state=c.state,kwargs...)=
            calculate_ace(props,sys,c;parameters,state,kwargs...)
        AC.calculate(props::$PropertyType,sys::AbstractSystem{3},c::ACECalculator,parameters,state;kwargs...)=
            calculate_ace(props,sys,c;parameters,state,kwargs...)
    end
end
