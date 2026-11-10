using PotentialFrameWork, Test, Unitful, AtomsBase, SpheriCart, LinearAlgebra
include(joinpath(ENV["PFW_ROOT"],"test/fixtures.jl"))
const ACECore=PotentialFrameWork.ACE
const ACE_MODEL_PATH=joinpath(ENV["PFW_ROOT"],"data/ace_ref/resolved_model_v1.toml")
include(joinpath(ENV["PFW_ROOT"],"test/ace_products.jl"))
