pushfirst!(LOAD_PATH,normpath(joinpath(@__DIR__,"../..")))
using PotentialFrameWork, InteractiveUtils
const A=PotentialFrameWork.ACE
model,ps,st=A.load_resolved(joinpath(@__DIR__,"../../data/ace_ref/resolved_model_v1.toml"))
R=zeros(4,length(model.radial_spec)); dR=similar(R)
args=(R,dR,(1.1,1.7,2.3,3.1),model.radial[1,2],1,Val(4),(true,true,true,true))
open(joinpath(@__DIR__,"radial-llvm-julia$(VERSION).txt"),"w") do io
    code_llvm(io,A.radial_block!,typeof(args);optimize=true)
end
println("Recorded optimized LLVM for width-4 joint radial values/derivatives")
