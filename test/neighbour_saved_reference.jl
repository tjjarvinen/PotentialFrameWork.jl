# Standalone saved-oracle regression; no ACEpotentials or reference regeneration.
# julia --project=<application-env> --threads=4 test/neighbour_saved_reference.jl
using PotentialFrameWork, AtomsBase, Unitful, SpheriCart, TOML, SHA, Test
include("fixtures.jl")
const A=PotentialFrameWork.ACE
fixture_dir=joinpath(@__DIR__,"../data/ace_ref")
reference=TOML.parsefile(joinpath(fixture_dir,"reloaded_predictions.toml"))
model,ps,st=A.load_resolved(joinpath(fixture_dir,"resolved_model_v1.toml"))
frames=TestSystem[]
# Deliberately small reader for this fixture's declared species/pos columns and
# periodic Lattice header. No topology comes from the saved oracle.
open(joinpath(fixture_dir,"reference_results.xyz")) do io
    while !eof(io)
        n=parse(Int,readline(io)); header=readline(io)
        lattice=match(r"Lattice=\"([^\"]+)\"",header)
        isnothing(lattice) && error("fixture lacks lattice")
        occursin("pbc=[T, T, T]",header) || error("reader supports this fully periodic fixture only")
        H=reshape(parse.(Float64,split(lattice[1])),3,3)
        X=zeros(3,n); labels=Symbol[]
        for i=1:n
            fields=split(readline(io)); push!(labels,Symbol(fields[1]))
            X[:,i]=parse.(Float64,fields[2:4])
        end
        push!(frames,TestSystem(X,H,(true,true,true),labels))
    end
end
@testset "neighborhood CPU vs saved full-precision 33-frame oracle" begin
    @test bytes2hex(sha256(read(joinpath(fixture_dir,"reference_model.json"))))==reference["model_sha256"]
    @test bytes2hex(sha256(read(joinpath(fixture_dir,"reference_results.xyz"))))==reference["geometry_sha256"]
    @test length(frames)==length(reference["frames"])==33
    for width in (1,4),batch_size in (1,2,5)
        env=A.FullEnvironment(;width); executor=A.CPUExecutor(model,env;slots=4,batch_size)
        for (i,sys) in enumerate(frames)
            input=A.ace_inputs(sys,model;environment=env)
            A.prepare_executor!(executor,env)
            A.execute_prepared!(executor,input.positions,input.cell,ps;parallel=true)
            y=A.prepared_result(executor); ref=reference["frames"][i]
            @test y.energy ≈ ref["energy_eV"] atol=1e-8 rtol=0
            @test y.forces ≈ hcat(ref["forces_eV_per_A"]...) atol=1e-9 rtol=1e-10
            @test y.virial ≈ hcat(ref["virial_eV_columns"]...) atol=1e-8 rtol=1e-10
        end
    end
end
