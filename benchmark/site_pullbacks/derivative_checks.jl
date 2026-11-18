using PotentialFrameWork, SpheriCart, ChainRulesCore, ForwardDiff, LinearAlgebra
include("../../examples/site_training_support.jl")
let
    A=PotentialFrameWork.ACE; NL=PotentialFrameWork.NeighbourList; PP=PotentialFrameWork.PairPotentials
    l,fp,n=SiteTrainingSupport.setup(;periodic=true)
    ps=(features=fp,head=(weights=[0.7,-0.3,0.4,0.6],bias=0.2,amplitude=0.8,quadratic=0.11))
    model=A.ACEFeatureEnergy(l,SiteTrainingSupport.SmoothInvariantHead()); X=n.positions
    bars=(energy=0.3,forces=reshape(cos.(1:9),3,3),virial=reshape(sin.(1:9),3,3),features=nothing)
    function loss(x,p)
        r=A.feature_energy_evaluate(x,model,p,n)
        bars.energy*r.energy+sum(bars.forces.*r.forces)+sum(bars.virial.*r.virial)
    end
    g=NL.site_pullback(X,A.FeatureEnergySiteModel(model),ps,n,bars)
    flat=PP.parameter_flatten(ps); v=sin.(1:length(flat)); V=reshape(cos.(1:9),3,3)
    predicted_p=dot(PP.parameter_flatten(g.parameters),v); predicted_x=sum(g.positions.*V)
    for h in (1e-4,1e-5,2e-6)
        measured_p=(loss(X,PP.parameter_rebuild(ps,flat+h*v))-loss(X,PP.parameter_rebuild(ps,flat-h*v)))/(2h)
        measured_x=(loss(X+h*V,ps)-loss(X-h*V,ps))/(2h)
        println("h=",h," parameter abs error=",abs(predicted_p-measured_p)," position abs error=",abs(predicted_x-measured_x))
    end
    blocks=A.ace_features(X,l,fp,n.environment,n.cell)
    hb=map(x -> reshape(sin.(1:length(x)),size(x)),blocks)
    jv=A.feature_jvp(X,l,fp,n.environment,n.cell,V)
    vjp=A.feature_vjp(X,l,fp,n.environment,n.cell,hb)
    println("feature JVP/VJP adjoint abs error=",abs(sum(sum(a.*b) for (a,b) in zip(hb,jv))-sum(vjp.positions.*V)))
end
