# julia --project=<test-env> -t4 benchmark/site_pullbacks/diagnostics.jl
using PotentialFrameWork, AtomsBase, StaticArrays, SpheriCart, ChainRulesCore, ForwardDiff, LinearAlgebra
include("../../examples/site_training_support.jl")
const A=PotentialFrameWork.ACE
const NL=PotentialFrameWork.NeighbourList
const NI=PotentialFrameWork.NeighbourInterface
function measurement(name,f)
    f(); f()
    r=@timed f()
    println(name,": seconds=",r.time,", bytes=",r.bytes)
    r.value
end
let
    l,fp,n=SiteTrainingSupport.setup()
    model=A.ACEFeatureEnergy(l,SiteTrainingSupport.SmoothInvariantHead())
    ps=(features=fp,head=(weights=[0.7,-0.3,0.4,0.6],bias=0.2,amplitude=0.8,quadratic=0.11))
    p=A.FeatureEnergySiteModel(model); X=n.positions
    b=(energy=0.3,forces=reshape(sin.(1:9),3,3),virial=reshape(cos.(1:9),3,3),features=nothing)
    measurement("generated basis + geometry preparation",()->SiteTrainingSupport.setup())
    sys=SiteTrainingSupport.TrainingSystem(copy(X),Matrix(4.0I,3,3),(false,false,false),ChemicalSpecies.([:Si,:O,:Si]))
    measurement("geometry preparation only",()->NL.prepare_neighbors(sys,l))
    measurement("composed chunk forward",()->NL.site_evaluate(X,p,ps,n;scatter_capacity=3,batch_size=3))
    y,pb=measurement("forward + pullback snapshot",()->rrule(NL.site_evaluate,X,p,ps,n;scatter_capacity=3,batch_size=3))
    measurement("combined backward + owning export",()->NL.site_pullback(X,p,ps,n,b;scatter_capacity=3))
    measurement("parameter-only backward + owning export",()->NL.site_pullback(X,p,ps,n,b;scatter_capacity=3,request=NI.SiteInputs(false,true)))
    measurement("backward from saved snapshot",()->pb(b))
    owned=deepcopy((X,p,ps,n))
    println("saved mathematical inputs bytes=",Base.summarysize(owned))
    capacity=maximum((n.environment.offsets[end,i]-n.environment.offsets[1,i] for i in axes(X,2));init=0)
    println("reusable parameter reverse trees bytes=",Base.summarysize(NI.make_scratch(p,capacity,Float64,ps)))
    println("global parameter accumulator bytes=",Base.summarysize(NI.parameter_buffer(p,ps,Float64)))
    println("atom positions gradient bytes=",sizeof(X),", scatter bytes=",3*sizeof(SVector{3,Float64}))
    # Separate export cost from prepared executor execution.
    w=NL.form_numa_workers(n,1;scatter_capacity=3)
    NL.foreach_worker!(w,1;scatter=(x=SVector{3,Float64},)) do worker end
    measurement("original-order export only",()->NL.get_results(w))
    m,sp,st=A.load_resolved(joinpath(@__DIR__,"../../data/ace_ref/resolved_model_v1.toml"))
    sp=A.trainable_parameters(m,sp)
    # Reuse numerical geometry with explicit Ti/Al channel identities.
    sn=deepcopy(n)
    scalar=A.ScalarSiteModel(m)
    measurement("scalar functional forward",()->A.ace_evaluate(X,m,sp,sn.environment,sn.cell))
    measurement("scalar chunk forward",()->NL.site_evaluate(X,scalar,sp,sn))
    measurement("scalar functional combined VJP",()->A.ace_vjp(X,m,sp,sn.environment,sn.cell,b.energy,b.forces,b.virial))
    measurement("scalar chunk combined VJP",()->NL.site_pullback(X,scalar,sp,sn,b))
    println("Julia=",VERSION,", threads=",Threads.nthreads(),", CPU=",Sys.CPU_NAME)
end
