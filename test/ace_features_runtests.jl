# Standalone entrypoint; the regular suite also includes ace_features.jl.
using PotentialFrameWork, Test, Unitful, AtomsBase
include("fixtures.jl")
const PP=PotentialFrameWork.PairPotentials
const ACECore=PotentialFrameWork.ACE
const ACE_MODEL_PATH=joinpath(@__DIR__,"../data/ace_ref/resolved_model_v1.toml")
include("ace_features.jl")
println("Feature validation: Julia ",VERSION,", threads=",Threads.nthreads())
