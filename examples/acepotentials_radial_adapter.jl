# Run only in the unchanged pinned environment:
# julia --project=data/ace_ref examples/acepotentials_radial_adapter.jl
using ACEpotentials
pushfirst!(LOAD_PATH,normpath(joinpath(@__DIR__,"..")))
using PotentialFrameWork, ExtXYZ, Test
const A=PotentialFrameWork.ACE
VERSION==v"1.13.2" && pkgversion(ACEpotentials)==v"0.10.2" || error("use pinned reference versions")
struct PinnedRadial{B}
    basis::B
    central::Int
    neighbour::Int
end
function (f::PinnedRadial)(v,d,r,ps,config)
    values,derivatives=ACEpotentials.Models.evaluate_ed(f.basis,r,f.central,f.neighbour,NamedTuple(),NamedTuple())
    copyto!(v,values); isnothing(d) || copyto!(d,derivatives)
    nothing
end
oracle,_=ACEpotentials.load_model(joinpath(@__DIR__,"../data/ace_ref/reference_model.json"))
model,ps,st=A.load_resolved(joinpath(@__DIR__,"../data/ace_ref/resolved_model_v1.toml"))
# evaluate_ed includes the original transform, spline normalization AND envelope.
# This allocating scalar adapter claims neither SIMD nor training AD support,
# and has no portable persistence registration. Radial compatibility is narrower
# than compatibility of arbitrary full ACEpotentials models.
metadata=A.RadialMetadata(model.radial_spec;cutoff=5.5,index_weights=sum.(model.radial_spec))
channels=[A.CallableRadial(PinnedRadial(oracle.model.rbasis,zi,zj),metadata) for zi in model.atomic_numbers,zj in model.atomic_numbers]
@testset "pinned ACEpotentials scalar radial adapter" begin
    for i=1:2,j=1:2,r in (0.31,1.17,2.73,5.21)
        a,da=A.radial_ed(channels[i,j],r)
        b,db=A.radial_ed(model.radial[i,j],r)
        @test a ≈ b atol=1e-11 rtol=1e-11
        @test da ≈ db atol=1e-10 rtol=1e-10
    end
    adapted=A.replace_radial(model,ps,channels)
    sys=first(ExtXYZ.load(joinpath(@__DIR__,"../data/ace_ref/reference_results.xyz")))
    input=A.ace_inputs(sys,model)
    expected=A.ace_evaluate(input.positions,model,ps,input.environment,input.cell)
    actual=A.ace_evaluate(input.positions,adapted.model,adapted.parameters,input.environment,input.cell)
    @test actual.energy ≈ expected.energy atol=1e-9
    @test actual.forces ≈ expected.forces atol=1e-9
    @test actual.virial ≈ expected.virial atol=1e-9
end
