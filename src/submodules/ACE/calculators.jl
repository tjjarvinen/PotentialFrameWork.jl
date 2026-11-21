"""Convenient owning inference wrapper. `parameters` is authoritative; update it
with `AtomsCalculators.set_parameters!`. Environment is rebuilt per public call.
Convenience calls have one active call per calculator; concurrent explicit calls
use distinct ACEEvaluationState instances.
"""
mutable struct ACECalculator{M,P,S,W}
    model::M
    parameters::P
    state::S
    environment::FullEnvironment
    executor::W
    evaluation_state::Any
end
ACECalculator(m,p,s,e,w)=ACECalculator(m,p,s,e,w,nothing)
ACECalculator{M,P,S,W}(m,p,s,e,w) where {M,P,S,W}=ACECalculator{M,P,S,W}(m,p,s,e,w,nothing)
function ACECalculator(model,ps,st=NamedTuple();width=4,slots=Threads.nthreads(),precision=Float64,angular_backend=:sphericart,batch_size=1,array_type=Array,product_algorithm=:auto)
    if array_type!==Array || model isa GPUACEModel
        b=adapt_model(model,ps;array_type,precision)
        if b.model isa GPUACEModel
            return GPUACECalculator(b.model,b.parameters,deepcopy(st),nothing)
        end
        model=b.model; ps=b.parameters
    elseif (eltype(model.reference)!==precision || !coherent_float_type(ps,precision)) && parameter_seed(ps) isa Union{Float32,Float64,Integer}
        b=adapt_model(model,ps;precision); model=b.model; ps=b.parameters
    end
    env=FullEnvironment(;width)
    executor=CPUExecutor(model,env;slots,precision,angular_backend,batch_size,product_algorithm)
    ACECalculator(model,deepcopy(ps),deepcopy(st),env,executor)
end
"""Join the calculator's default chunk workers. Explicit states close separately."""
function Base.close(c::ACECalculator)
    lock(ACE_STATE_INIT_LOCK) do
        isnothing(c.evaluation_state) || close(c.evaluation_state)
    end
    nothing
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
function calculator_result(sys,c;parameters=c.parameters,state=nothing,energy_only=false,parallel=false,executor=(parallel && !(state isa ACEEvaluationState) ? :colored : :chunks),simd=true,fused=false,product_algorithm=c.executor.prepared.product_algorithm,centers=nothing,properties=(AC.Energy(),AC.Forces(),AC.Virial()),
                           workers=Threads.nthreads(:default),groups=1,reuse_topology=false,timings=nothing,kwargs...)
    parameters=(isnothing(parameters) || parameters===NamedTuple()) ? c.parameters : parameters
    executor in (:chunks,:threaded,:colored) || throw(ArgumentError("executor must be :chunks, :threaded or :colored"))
    if executor===:threaded || (executor===:chunks && (!energy_only || state isa ACEEvaluationState))
        groups isa Integer && groups>0 || throw(ArgumentError("groups must be a positive integer"))
        workers isa Integer && workers>0 || throw(ArgumentError("workers must be a positive integer"))
        explicit=state isa ACEEvaluationState
        shared=isnothing(state) || (executor===:threaded && state===c.state)
        evaluation_state=explicit ? state : shared ? default_evaluation_state!(c) : ACEEvaluationState(c;state)
        try
            return chunk_calculator_result(sys,c,evaluation_state,parameters;energy_only,parallel,simd,fused,
                product_algorithm,centers,properties,workers=executor===:threaded ? workers : 1,
                groups,reuse_topology,timings,kwargs...)
        finally
            explicit || shared || close(evaluation_state)
        end
    end
    reuse_topology && throw(ArgumentError("topology reuse requires a chunk executor and prepared evaluation state"))
    state isa ACEEvaluationState && throw(ArgumentError("explicit ACEEvaluationState requires a chunk executor"))
    fast_parameters=parameter_seed(parameters) isa Union{Float32,Float64}
    boundary=fast_parameters && eltype(c.executor.force)===Float32 ? precision_geometry(sys,Float32) : sys
    input=fast_parameters ? ace_inputs(boundary,c.model;environment=c.environment,precision=eltype(c.executor.force)) : ace_inputs(sys,c.model;environment=c.environment)
    y=if energy_only || !fast_parameters
        ace_evaluate(input.positions,c.model,parameters,input.environment,input.cell;state=isnothing(state) ? c.state : state,energy_only,centers)
    else
        prepare_executor!(c.executor,c.environment)
        execute_prepared!(c.executor,input.positions,input.cell,parameters;parallel,simd,fused,product_algorithm,centers,properties)
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
    r=calculator_result(sys,c;parameters,state,energy_only=all(p -> p isa AC.Energy,props),properties=normalized,kwargs...)
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
