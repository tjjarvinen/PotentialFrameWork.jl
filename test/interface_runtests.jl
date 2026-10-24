# Focused reproduction in an application environment containing the test weakdeps:
# julia --project=<env> --threads=1 test/interface_runtests.jl
using PotentialFrameWork, Test, Unitful, AtomsBase, SpheriCart, LinearAlgebra
using ForwardDiff, Zygote, ChainRulesCore, LuxCore
import AtomsCalculators as AC
const ACECore=PotentialFrameWork.ACE
const PP=PotentialFrameWork.PairPotentials
const ACE_MODEL_PATH=joinpath(@__DIR__,"../data/ace_ref/resolved_model_v1.toml")
include("fixtures.jl")
function central_gradient(f,x;h=2e-6)
    map(eachindex(x)) do k
        xp=copy(x); xm=copy(x); xp[k]+=h; xm[k]-=h
        (f(xp)-f(xm))/(2h)
    end |> g -> reshape(g,size(x))
end
include("ace_interface.jl")
println("Interface validation: Julia ",VERSION,", threads=",Threads.nthreads())
