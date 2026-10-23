using ACEpotentials
pushfirst!(LOAD_PATH,normpath(joinpath(@__DIR__,"../..")))
using PotentialFrameWork, ACEpotentials, ExtXYZ, Test, LinearAlgebra
const A=PotentialFrameWork.ACE
const PP=PotentialFrameWork.PairPotentials
const Z=ACEpotentials.Models.Zygote
model,readout,st=A.load_resolved(joinpath(@__DIR__,"resolved_model_v1.toml"))
ps=A.trainable_parameters(model,readout)
sys=first(ExtXYZ.load(joinpath(@__DIR__,"reference_results.xyz")))
inp=A.ace_inputs(sys,model;environment=A.FullEnvironment(width=4))
X=inp.positions; H=inp.cell; env=inp.environment
B=sin.(X); C=reshape(cos.(1:9),3,3)
@testset "periodic multi-image ACE derivative products" begin
    for kind in (:energy,:forces,:virial,:combined)
        function loss(x,p)
            r=A.ace_evaluate(x,model,p,env,H)
            kind===:energy && return r.energy
            kind===:forces && return sum(B.*r.forces)
            kind===:virial && return sum(C.*r.virial)
            0.7r.energy+sum(B.*r.forces)+sum(C.*r.virial)
        end
        gx,gp=Z.gradient(loss,X,ps)
        h=2e-6
        fd=map(eachindex(X)) do k
            plus=copy(X); minus=copy(X); plus[k]+=h; minus[k]-=h
            (loss(plus,ps)-loss(minus,ps))/(2h)
        end
        @test vec(gx) ≈ fd atol=3e-6 rtol=3e-6
        flat=PP.parameter_flatten(ps)
        for key in keys(ps)
            direction=PP.parameter_rebuild(ps,zeros(length(flat)))
            target=direction[key] isa NamedTuple ? direction[key].scale : direction[key]
            target .= reshape(sin.(1:length(target)),size(target))
            d=PP.parameter_flatten(direction)
            fd=(loss(X,PP.parameter_rebuild(ps,flat+h*d))-loss(X,PP.parameter_rebuild(ps,flat-h*d)))/(2h)
            @test dot(PP.parameter_flatten(gp),d) ≈ fd atol=3e-6 rtol=3e-6
        end
    end
end
