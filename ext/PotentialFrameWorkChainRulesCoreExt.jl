module PotentialFrameWorkChainRulesCoreExt
using PotentialFrameWork, ChainRulesCore
import ChainRulesCore: rrule

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
            project(map(x -> x isa Real ? delta*x : delta.*x,g))
        end
        (NoTangent(),@not_implemented("system derivatives are outside parameter-only energy AD"),
         @not_implemented("differentiate explicit parameters, not the bundled potential"),dp)
    end
    y,energy_pullback
end
end
