# julia --project=<test-env> examples/ace_feature_training.jl
using PotentialFrameWork, SpheriCart, ChainRulesCore, ForwardDiff, Zygote, LuxCore
include("site_training_support.jl")
using .SiteTrainingSupport: setup, step
let
    A=PotentialFrameWork.ACE
    layer,ps,n=setup()
    lux=lux_layer(layer;parameters=ps,executor=:chunks)
    input=(positions=n.positions,neighbors=n,batch_size=3,scatter_capacity=3)
    teacher=deepcopy(ps)
    for W in values(teacher.mixing); W .*= 1.4; end
    target=first(LuxCore.apply(lux,input,teacher,NamedTuple()))
    loss(p)=sum(sum(abs2,a-b) for (a,b) in zip(first(LuxCore.apply(lux,input,p,NamedTuple())),target))
    before=loss(ps)
    for _=1:40
        g=Zygote.gradient(loss,ps)[1]
        @assert all(isfinite,PotentialFrameWork.PairPotentials.parameter_flatten(g))
        ps=step(ps,g,0.1)
    end
    after=loss(ps)
    @assert after<before
    println("Direct feature training: steps=40, eta=0.1, loss=",before," -> ",after)
end
