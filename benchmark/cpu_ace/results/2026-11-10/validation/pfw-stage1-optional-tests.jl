using PotentialFrameWork, Test, SpheriCart, ACEfit, LinearAlgebra, Unitful, AtomsBase
const ACE_MODEL_PATH=joinpath(ENV["PFW_ROOT"],"data/ace_ref/resolved_model_v1.toml")
source=read(joinpath(ENV["PFW_ROOT"],"test/ace_interface_optional.jl"),String)
include_string(Main,replace(source,"include(\"ace_interface_plotting.jl\")"=>""),"ace_interface_optional_without_plotting.jl")
