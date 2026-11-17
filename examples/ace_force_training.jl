# julia --project=<test-env> examples/ace_force_training.jl
using PotentialFrameWork, SpheriCart, ChainRulesCore, ForwardDiff, Zygote, LuxCore
include("site_training_support.jl")
using .SiteTrainingSupport: setup, step, SmoothInvariantHead
let
    A=PotentialFrameWork.ACE; PP=PotentialFrameWork.PairPotentials
    layer,fp,n=setup()
    model=A.ACEFeatureEnergy(layer,SmoothInvariantHead())
    ps=(features=fp,head=(weights=[0.7,-0.3,0.4,0.6],bias=0.2,amplitude=0.8,quadratic=0.11))
    teacher=deepcopy(ps)
    for W in values(teacher.features.mixing); W .*= 1.3; end
    teacher.head.weights .*= 1.2
    lux=lux_layer(model;parameters=ps)
    input=(positions=n.positions,neighbors=n,batch_size=3,scatter_capacity=3)
    target=first(LuxCore.apply(lux,input,teacher,NamedTuple()))
    function loss(p)
        r=first(LuxCore.apply(lux,input,p,NamedTuple()))
        (r.energy-target.energy)^2+10sum(abs2,r.forces-target.forces)
    end
    before=loss(ps)
    before_prediction=first(LuxCore.apply(lux,input,ps,NamedTuple()))
    for _=1:60
        g=Zygote.gradient(loss,ps)[1]
        @assert all(isfinite,PP.parameter_flatten(g))
        ps=step(ps,g,0.03)
    end
    after=loss(ps)
    @assert after<before
    after_prediction=first(LuxCore.apply(lux,input,ps,NamedTuple()))
    println("Force MSE sum: ",sum(abs2,before_prediction.forces-target.forces)," -> ",sum(abs2,after_prediction.forces-target.forces))
    println("Composed E/F training: steps=60, eta=0.03, force weight=10, loss=",before," -> ",after)
end
