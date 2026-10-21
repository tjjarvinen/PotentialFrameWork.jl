module PotentialFrameWorkForwardDiffExt
using PotentialFrameWork, ForwardDiff
import PotentialFrameWork: energy_parameter_gradient

# The first parameter interface supports flat NamedTuples of real scalars/vectors.
# Preserve field names and shapes at the ChainRules boundary.
function unpack(template::NamedTuple, flat)
    offset=0
    vals=map(values(template)) do x
        if x isa Real
            offset+=1
            flat[offset]
        else
            n=length(x)
            v=flat[offset+1:offset+n]
            offset+=n
            v
        end
    end
    NamedTuple{keys(template)}(vals)
end
function energy_parameter_gradient(sys,p::AbstractPairPotential,parameters::NamedTuple;kwargs...)
    PotentialFrameWork.PairPotentials.parameter_seed(parameters)
    flat=vcat(map(x -> x isa Real ? [x] : collect(x),values(parameters))...)
    isempty(flat) && return unpack(parameters,flat)
    gradient=ForwardDiff.gradient(flat) do theta
        parameter_energy(sys,p,unpack(parameters,theta);kwargs...)
    end
    unpack(parameters,gradient)
end
end
