# Serial correctness-path costs. Run in an application environment containing
# SpheriCart, ChainRulesCore and ForwardDiff, with this checkout developed.
# julia --project=<env> -t4 benchmark/ace_features.jl
using PotentialFrameWork, SpheriCart, ChainRulesCore, ForwardDiff, AtomsBase, Unitful
using LinearAlgebra
const A=PotentialFrameWork.ACE

function feature_costs()
    radial=A.CallableRadial(A.SmoothExponential(),
        A.RadialMetadata([(1,0),(1,1),(2,1)];cutoff=4.,ad=true);
        config=(cutoff=4.,powers=[0,0,1]))
    layer=A.ACEFeatureLayer(species=[:Si,:O],radial=radial,correlation=2,lmax=1,
        output_irreps=((L=0,parity=:even,channels=3),
                      (L=1,parity=:odd,channels=2),(L=1,parity=:even,channels=2)))
    ps=A.feature_parameters(layer;custom_radial=(rate=[0.4],amplitude=[1.3]))
    system=AtomsBase.isolated_system([:Si=>[0.1,0.2,0.3]*u"Å",
        :O=>[1.2,0.4,1.1]*u"Å",:Si=>[0.7,1.8,0.5]*u"Å"])
    input=A.ace_inputs(system,layer)
    forward()=A.ace_features(input.positions,layer,ps,input.environment,input.cell)
    rule()=rrule(A.ace_features,input.positions,layer,ps,input.environment,input.cell)
    y,pb=rule(); cot=map(b -> ones(size(b)),y)
    forward(); rule(); pb(cot)
    GC.gc()
    println("Julia=",VERSION," CPU=",Sys.CPU_NAME," threads=",Threads.nthreads(),
        " BLAS=",BLAS.get_num_threads())
    println("basis=",A.feature_basis_summary(layer)," density=",length(layer.density),
        " products=",length(layer.products))
    println("forward_allocated_bytes=",@allocated forward())
    println("rrule_forward_allocated_bytes=",@allocated rule())
    println("pullback_allocated_bytes=",@allocated pb(cot))
    println("output_payload_bytes=",sum(sizeof,y),
        " saved_pullback_summarysize=",Base.summarysize(pb))
    println("median_forward_seconds=",sort([@elapsed forward() for _=1:7])[4])
    println("median_pullback_seconds=",sort([@elapsed pb(cot) for _=1:7])[4])
    for m in (SpheriCart,ChainRulesCore,ForwardDiff)
        println(nameof(m),"=",pkgversion(m))
    end
end
feature_costs()
