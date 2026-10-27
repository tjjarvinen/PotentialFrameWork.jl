module PotentialFrameWorkChainRulesCoreExt
using PotentialFrameWork, ChainRulesCore
import ChainRulesCore: rrule
using PotentialFrameWork.ACE: ace_evaluate, ace_vjp, ResolvedModel, FullEnvironment

function rrule(::typeof(parameter_energy),sys,p::AbstractPairPotential,parameters::NamedTuple;kwargs...)
    y=parameter_energy(sys,p,parameters;kwargs...)
    # Capture derivatives now: a supplied mutable workspace, system or parameter
    # vector may be reused/changed before the pullback is called.
    g=energy_parameter_gradient(sys,p,parameters;kwargs...)
    project=ProjectTo(parameters)
    function energy_pullback(delta)
        delta=unthunk(delta)
        dp=if delta isa AbstractZero
            ZeroTangent()
        else
            project(PotentialFrameWork.PairPotentials.parameter_scale(delta,g))
        end
        (NoTangent(),@not_implemented("system derivatives are outside parameter-only energy AD"),
         @not_implemented("differentiate explicit parameters, not the bundled potential"),dp)
    end
    y,energy_pullback
end

function rrule(::typeof(pair_evaluate),X,p::AbstractPairPotential,ps::NamedTuple,
               env::PairEnvironment,H;request=EnergyForcesVirial(),state=NamedTuple())
    y=pair_evaluate(X,p,ps,env,H;request,state)
    # Own every numerical input and discrete map needed for recomputation. Neither
    # caller mutation nor reuse of an unrelated inference workspace changes this.
    saved=deepcopy((X,p,ps,env,H))
    projectX=ProjectTo(X); projectps=ProjectTo(ps)
    function pair_pullback(delta)
        delta=unthunk(delta)
        if delta isa AbstractZero
            return (NoTangent(),ZeroTangent(),NoTangent(),ZeroTangent(),NoTangent(),NoTangent())
        end
        cot(name,zero_value)=if !hasproperty(delta,name)
            zero_value
        else
            v=unthunk(getproperty(delta,name))
            v isa AbstractZero ? zero_value : v
        end
        de=cot(:energy,zero(y.energy))
        df=cot(:forces,zeros(eltype(X),size(X)))
        dw=cot(:virial,zeros(eltype(X),3,3))
        g=pair_vjp(saved...,de,df,dw)
        (NoTangent(),projectX(g.positions),NoTangent(),projectps(g.parameters),NoTangent(),
         @not_implemented("cell derivatives are outside the fixed-cell pair AD contract"))
    end
    y,pair_pullback
end

function rrule(::typeof(ace_evaluate),X,model::ResolvedModel,ps,env::FullEnvironment,H;
               state=NamedTuple(),energy_only=false,centers=nothing)
    y=ace_evaluate(X,model,ps,env,H;state,energy_only,centers)
    saved_centers=deepcopy(centers)
    saved=deepcopy((X,model,ps,env,H))
    px=ProjectTo(X); pp=ProjectTo(ps)
    function ace_pullback(delta)
        delta=unthunk(delta)
        if delta isa AbstractZero
            return (NoTangent(),ZeroTangent(),NoTangent(),ZeroTangent(),NoTangent(),NoTangent())
        end
        function cot(name,z)
            hasproperty(delta,name) || return z
            v=unthunk(getproperty(delta,name))
            v isa AbstractZero ? z : v
        end
        de=cot(:energy,zero(y.energy)); df=cot(:forces,zeros(eltype(X),size(X))); dw=cot(:virial,zeros(eltype(X),3,3))
        g=ace_vjp(saved...,de,df,dw;centers=saved_centers)
        (NoTangent(),px(g.positions),NoTangent(),pp(g.parameters),NoTangent(),
         @not_implemented("cell derivatives are outside the fixed-cell ACE AD contract"))
    end
    y,ace_pullback
end
# These entrypoints deliberately implement inference only.
for f in (PotentialFrameWork.PairPotentials.evaluate_gpu,
          PotentialFrameWork.PairPotentials.evaluate_gpu!,
          PotentialFrameWork.PairPotentials.execute_gpu!,
          PotentialFrameWork.PairPotentials.prepare_gpu_execution,
          PotentialFrameWork.NeighbourList.gpu_execute!)
    @eval function ChainRulesCore.rrule(::typeof($f),args...;kwargs...)
        throw(ArgumentError("GPU differentiation is unsupported; use the functional CPU pair_evaluate path"))
    end
    @eval function ChainRulesCore.frule(tangents,::typeof($f),args...;kwargs...)
        throw(ArgumentError("GPU differentiation is unsupported; use the functional CPU pair_evaluate path"))
    end
end

end
