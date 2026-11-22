module SitePullbacksTests
using PotentialFrameWork, Test, Unitful, LuxCore, StaticArrays, LinearAlgebra, ForwardDiff, ChainRulesCore, Zygote, SpheriCart
include("fixtures.jl")
import AtomsCalculators as AC
const NI=PotentialFrameWork.NeighbourInterface
const NL=PotentialFrameWork.NeighbourList
const A=PotentialFrameWork.ACE
const PP=PotentialFrameWork.PairPotentials
struct NonlinearSite end
NI.make_scratch(::NonlinearSite,c,T)=nothing
mutable struct NonlinearOutput{T}
    energy::T
    gradients::Vector{SVector{3,T}}
    features::Tuple{Matrix{T}}
end
NI.make_site_output(::NonlinearSite,c,T,o)=NonlinearOutput(zero(T),zeros(SVector{3,T},c),(zeros(T,2,1),))
NI.output_energy(o::NonlinearOutput)=o.energy
NI.output_gradient(o::NonlinearOutput,k)=o.gradients[k]
NI.output_features(o::NonlinearOutput)=o.features
NI.site_derivative_capabilities(::NonlinearSite)=(energy=NI.SiteInputs(),gradients=NI.SiteInputs(),features=NI.SiteInputs())
function NI.evaluate_site!(o,p::NonlinearSite,s,ps,scratch,request::NI.SiteOutputs)
    q=sum((sum(abs2,NI.displacement(s,k)) for k=1:NI.neighbour_count(s) if NI.active(s,k));init=zero(o.energy))
    a=ps.a; b=ps.nested.b[1]
    o.energy=a^2*q/2+b*q^2/4
    for k=1:NI.neighbour_count(s)
        NI.active(s,k) && (o.gradients[k]=(a^2+b*q)*NI.displacement(s,k))
    end
    o.features[1][:,1]=[q,a*q^2]
    nothing
end
function NI.site_vjp!(rb,pb,p::NonlinearSite,s,ps,bars,scratch;request=NI.SiteInputs())
    edges=[k for k=1:NI.neighbour_count(s) if NI.active(s,k)]
    q=sum((sum(abs2,NI.displacement(s,k)) for k in edges);init=0.)
    a=ps.a; b=ps.nested.b[1]
    de=isnothing(bars.energy) ? 0. : bars.energy
    c=isnothing(bars.gradients) ? 0. : sum((dot(bars.gradients[k],NI.displacement(s,k)) for k in edges);init=0.)
    h=isnothing(bars.features) ? zeros(2) : vec(bars.features[1])
    if request.displacements
        for k in edges
            r=NI.displacement(s,k)
            rb[k]+=(de*(a^2+b*q)+2b*c+2h[1]+4a*q*h[2])*r
            isnothing(bars.gradients) || (rb[k]+=(a^2+b*q)*bars.gradients[k])
        end
    end
    request.parameters && NI.add_parameter_buffer!(p,pb,(a=de*a*q+2a*c+h[2]*q^2,nested=(b=[de*q^2/4+q*c],)))
    nothing
end
function fd(f,x;h=2e-6)
    g=zeros(length(x))
    for k in eachindex(x)
        plus=copy(x); minus=copy(x); plus[k]+=h; minus[k]-=h
        g[k]=(f(plus)-f(minus))/(2h)
    end
    g
end

function contraction(r,b)
    (isnothing(b.energy) ? 0. : b.energy*r.energy)+(isnothing(b.forces) ? 0. : sum(b.forces.*r.forces))+
    (isnothing(b.virial) ? 0. : sum(b.virial.*r.virial))+
    (isnothing(b.features) ? 0. : sum(sum(sum(x.*y) for (x,y) in zip(u,v)) for (u,v) in zip(b.features,r.features)))
end
@testset "generic nonlinear chunk VJP" begin
    ps=(a=0.7,nested=(b=[0.13],)); model=NonlinearSite()
    sys=TestSystem([3.4 0.2 1.0 2.1; 0.1 0.5 1.2 0.2; 0.9 0.7 0.3 0.4],Matrix(5.0I,3,3),(false,false,false),fill(:H,4))
    for width in (1,4), centers in (nothing,[3,1],Int[])
        n=NL.prepare_neighbors(sys,2.2;width)
        @test n.p != collect(1:4)
        X=copy(n.positions); outputs=NI.SiteOutputs(true,true,true)
        y=NL.site_evaluate(X,model,ps,n;outputs,centers,scatter_capacity=4)
        bars=(energy=0.3,forces=reshape(sin.(1:12),3,4),virial=reshape(cos.(1:9),3,3),features=[(reshape([0.2i,-0.1i],2,1),) for i=1:4])
        g=NL.site_pullback(X,model,ps,n,bars;outputs,centers,scatter_capacity=4)
        @test vec(g.positions) ≈ fd(x -> contraction(NL.site_evaluate(reshape(x,3,4),model,ps,n;outputs,centers),bars),vec(X)) rtol=2e-6 atol=1e-7
        flat=PP.parameter_flatten(ps)
        @test PP.parameter_flatten(g.parameters) ≈ fd(t -> contraction(NL.site_evaluate(X,model,PP.parameter_rebuild(ps,t),n;outputs,centers),bars),flat) rtol=2e-6 atol=1e-7
        for request in (NI.SiteInputs(false,true),NI.SiteInputs(true,false))
            h=NL.site_pullback(X,model,ps,n,bars;outputs,centers,request)
            @test request.displacements ? h.positions≈g.positions : isnothing(h.positions)
            @test request.parameters ? PP.parameter_flatten(h.parameters)≈PP.parameter_flatten(g.parameters) : isnothing(h.parameters)
        end
        y,back=rrule(NL.site_evaluate,X,model,ps,n;outputs,centers)
        expected=back(bars)
        X .+= 0.1; ps.nested.b .+= 0.2
        @test back(bars)[2]≈expected[2]
        @test back(ZeroTangent())[2] isa AbstractZero
        ps.nested.b .-= 0.2
    end
end
@testset "scalar ACE chunk products" begin
    m,ps,st=A.load_resolved(joinpath(@__DIR__,"../data/ace_ref/resolved_model_v1.toml"))
    ps=A.trainable_parameters(m,ps)
    sys=TestSystem([0.1 1.6 2.3;0.2 0.4 1.1;0.3 0.7 0.2],Matrix(8.0I,3,3),(false,false,false),[:Ti,:Al,:Ti])
    calc=A.ACECalculator(m,ps,st); n=NL.prepare_neighbors(sys,calc;width=4)
    X=n.positions; model=A.ScalarSiteModel(m)
    for centers in (nothing,[2],Int[])
        r=NL.site_evaluate(X,model,ps,n;centers)
        ref=A.ace_evaluate(X,m,ps,n.environment,n.cell;centers)
        @test r.energy≈ref.energy
        @test r.forces≈ref.forces
        @test r.virial≈ref.virial
        for field in (:energy,:forces,:virial,:all)
            b=(energy=field in (:energy,:all) ? 0.3 : nothing,forces=field in (:forces,:all) ? reshape(sin.(1:9),3,3) : nothing,virial=field in (:virial,:all) ? reshape(cos.(1:9),3,3) : nothing,features=nothing)
            g=NL.site_pullback(X,model,ps,n,b;centers)
            refg=A.ace_vjp(X,m,ps,n.environment,n.cell,something(b.energy,0.),something(b.forces,zeros(3,3)),something(b.virial,zeros(3,3));centers)
            @test g.positions≈refg.positions atol=1e-8
            @test PP.parameter_flatten(g.parameters)≈PP.parameter_flatten(refg.parameters) atol=1e-8
        end
    end
end
@testset "chunk feature first derivatives" begin
    ch=A.CallableRadial(A.SmoothExponential(),A.RadialMetadata([(1,0),(1,1),(2,1),(1,2)];cutoff=4.,ad=true);config=(cutoff=4.,powers=[0,0,1,0]))
    l=A.ACEFeatureLayer(species=[:Si,:O],radial=ch,correlation=2,lmax=2,
        output_irreps=((L=0,parity=:even,channels=2),(L=1,parity=:odd,channels=2),(L=1,parity=:even,channels=1),(L=2,parity=:even,channels=2)))
    ps=A.feature_parameters(l;custom_radial=(rate=[0.4],amplitude=[1.3]))
    for (b,W) in enumerate(values(ps.mixing)); W .= reshape(sin.(b .+ (1:length(W))),size(W)); end
    sys=TestSystem([0.1 1.6 2.3;0.2 0.4 1.1;0.3 0.7 0.2],Matrix(8.0I,3,3),(false,false,false),[:Si,:O,:Si])
    n=NL.prepare_neighbors(sys,4.;channels=[1,2,1],species_count=2,width=4)
    X=n.positions; p=A.FeatureSiteModel(l); outputs=NI.SiteOutputs(false,false,true)
    y=NL.site_evaluate(X,p,ps,n;outputs,scatter_capacity=3)
    h=A.ace_features(X,l,ps,n.environment,n.cell)
    @test all(A.stack_feature_sites(y.features,l,3)[b]≈h[b] for b in eachindex(h))
    bars=map(x -> reshape(sin.(1:length(x)),size(x)),h)
    b=(energy=nothing,forces=nothing,virial=nothing,features=[map(x -> copy(x[:,:,i]),bars) for i=1:3])
    g=NL.site_pullback(X,p,ps,n,b;outputs)
    ref=A.feature_vjp(X,l,ps,n.environment,n.cell,bars)
    @test g.positions≈ref.positions
    @test PP.parameter_flatten(g.parameters)≈PP.parameter_flatten(ref.parameters)
    @test vec(g.positions)≈fd(x -> contraction(NL.site_evaluate(reshape(x,3,3),p,ps,n;outputs),b),vec(X)) rtol=3e-6 atol=1e-8
    flat=PP.parameter_flatten(ps); direction=sin.(1:length(flat)); eps=2e-6
    objective(t)=contraction(NL.site_evaluate(X,p,PP.parameter_rebuild(ps,t),n;outputs),b)
    @test dot(PP.parameter_flatten(g.parameters),direction)≈(objective(flat+eps*direction)-objective(flat-eps*direction))/(2eps) rtol=3e-6 atol=1e-8
    zg=Zygote.gradient(q -> sum(sum(abs2,block) for block in A.ace_features_chunks(X,l,q,n)),ps)[1]
    @test all(isfinite,PP.parameter_flatten(zg))
    _,back=rrule(A.ace_features_chunks,X,l,ps,n)
    for block in eachindex(h)
        cot=ntuple(k -> k==block ? bars[k] : zero(h[k]),length(h))
        truth=A.feature_vjp(X,l,ps,n.environment,n.cell,cot)
        lazy=ntuple(k -> k==block ? Thunk(() -> bars[k]) : ZeroTangent(),length(h))
        @test back(lazy)[2]≈truth.positions
    end
end
@testset "images empty buffers and batch fallback" begin
    model=NonlinearSite(); ps=(a=0.3,nested=(b=[0.01],)); outputs=NI.SiteOutputs(true,true,true)
    for count in (0,1,2,4), periodic in (false,true)
        sys=fixture(count;L=2.,mixed=false,pbc=(periodic,periodic,periodic))
        n=NL.prepare_neighbors(sys,2.1;width=4); X=n.positions
        r=NL.site_evaluate(X,model,ps,n;outputs,scatter_capacity=max(count,1),batch_size=3)
        @test r.energy≈NL.site_evaluate(X,model,ps,n;outputs,batch_size=1).energy
        b=(energy=0.2,forces=reshape(sin.(1:3count),3,count),virial=reshape(cos.(1:9),3,3),features=nothing)
        g=NL.site_pullback(X,model,ps,n,b;outputs,scatter_capacity=max(count,1))
        @test vec(g.positions)≈fd(x -> contraction(NL.site_evaluate(reshape(x,3,count),model,ps,n;outputs),b),vec(X)) rtol=3e-6 atol=1e-7
        @test sum(g.positions;dims=2)≈zeros(3,1) atol=1e-10
    end
    pbuf=NI.parameter_buffer(model,(empty=NamedTuple(),scalar=2.,array=zeros(2,3)),Float64)
    NI.add_parameter_buffer!(model,pbuf,(empty=NamedTuple(),scalar=3.,array=ones(2,3)))
    @test NI.parameter_value(model,pbuf).scalar==3.
    NI.reset_parameter_buffer!(model,pbuf)
    @test NI.parameter_value(model,pbuf).array==zeros(2,3)
    @test NI.parameter_value(model,NI.parameter_buffer(model,NamedTuple(),Float64))==NamedTuple()
end
struct EnergyOnlyDerivatives end
NI.site_derivative_capabilities(::EnergyOnlyDerivatives)=(energy=NI.SiteInputs(),gradients=NI.SiteInputs(false,false),features=NI.SiteInputs(false,false))
@testset "capability validation and thunk boundary" begin
    @test isnothing(NL.validate_site_derivatives(EnergyOnlyDerivatives(),NI.SiteOutputs(true,false,false),NI.SiteInputs()))
    @test_throws ArgumentError NL.validate_site_derivatives(EnergyOnlyDerivatives(),NI.SiteOutputs(false,true,false),NI.SiteInputs())
    sys=fixture(3;L=3.,mixed=false); n=NL.prepare_neighbors(sys,2.2;width=4)
    p=NonlinearSite(); ps=(a=0.3,nested=(b=[0.1],)); outputs=NI.SiteOutputs(true,true,true)
    y,pb=rrule(NL.site_evaluate,n.positions,p,ps,n;outputs,batch_size=3)
    plain=pb((energy=0.2,))
    thunk=pb(Thunk(() -> (energy=Thunk(() -> 0.2),forces=ZeroTangent(),features=ZeroTangent())))
    @test plain[2]≈thunk[2]
    @test plain[4].a≈thunk[4].a && plain[4].nested.b≈thunk[4].nested.b
    n.environment.shifts .+= Ref(SVector(1,0,0)); n.positions .+= 2.
    @test pb((energy=0.2,))[2]≈plain[2]
end
@testset "public chunk Lux and unitful calculator derivatives" begin
    m,ps,st=A.load_resolved(joinpath(@__DIR__,"../data/ace_ref/resolved_model_v1.toml"))
    ps=A.trainable_parameters(m,ps)
    sys=TestSystem([0.1 1.6;0.2 0.4;0.3 0.7],Matrix(8.0I,3,3),(false,false,false),[:Ti,:Al])
    c=A.ACECalculator(m,ps,st); n=NL.prepare_neighbors(sys,c)
    layer=PotentialFrameWork.lux_layer(m;parameters=ps,executor=:chunks)
    input=(positions=n.positions,neighbors=n,batch_size=3)
    g=Zygote.gradient(p -> sum(abs2,first(LuxCore.apply(layer,input,p,NamedTuple())).forces),ps)[1]
    @test all(isfinite,PP.parameter_flatten(g))
    state=A.ACEEvaluationState(c)
    objective(p)=ustrip(u"eV",AC.calculate(AC.Energy(),sys,c,p,state).energy)
    gp=Zygote.gradient(objective,ps)[1]
    reference=Zygote.gradient(p -> A.ace_evaluate_chunks(n.positions,m,p,n).energy,ps)[1]
    @test PP.parameter_flatten(gp)≈PP.parameter_flatten(reference)
    function unitloss(p)
        r=AC.calculate((AC.Energy(),AC.Forces(),AC.Virial()),sys,c,p,state)
        ustrip(u"eV",r.energy)*0.2+sum(sum(abs2,ustrip.(u"eV/Å",f)) for f in r.forces)+sum(abs2,ustrip.(u"eV",r.virial))*0.01
    end
    ug=Zygote.gradient(unitloss,ps)[1]
    ng=Zygote.gradient(ps) do p
        r=A.ace_evaluate_chunks(n.positions,m,p,n)
        0.2r.energy+sum(abs2,r.forces)+0.01sum(abs2,r.virial)
    end[1]
    @test PP.parameter_flatten(ug)≈PP.parameter_flatten(ng) rtol=1e-10
    fourarg=Zygote.gradient(p -> ustrip(u"eV",AC.calculate(AC.Energy(),sys,c,p).energy),ps)[1]
    @test PP.parameter_flatten(fourarg)≈PP.parameter_flatten(reference)
    threaded_sys=TestSystem([0.1 1.6 6.2;0.2 0.4 0.1;0.3 0.7 0.5],
                            Matrix(18.0I,3,3),(false,false,false),[:Ti,:Al,:Ti])
    thread_state=A.ACEEvaluationState(c)
    try
        function force_loss(p,executor)
            r=AC.calculate((AC.Energy(),AC.Forces(),AC.Virial()),threaded_sys,c,p,thread_state;executor,workers=2)
            ustrip(u"eV",r.energy)*0.2+sum(sum(abs2,ustrip.(u"eV/Å",f)) for f in r.forces)+sum(abs2,ustrip.(u"eV",r.virial))*0.01
        end
        _,back=Zygote.pullback(p -> force_loss(p,:threaded),ps)
        @test length(thread_state.chunks.pool.tasks)==2
        close(thread_state)
        serial=Zygote.gradient(p -> force_loss(p,:chunks),ps)[1]
        threaded=first(back(1.0))
        @test PP.parameter_flatten(threaded)≈PP.parameter_flatten(serial) rtol=1e-10
    finally
        close(thread_state)
    end
    state.active=true
    @test_throws ArgumentError Zygote.gradient(objective,ps)
    state.active=false; state.scatter_capacity=1
    @test_throws ArgumentError Zygote.gradient(objective,ps)
    @test !state.active
end
@testset "scalar projection callable and independent directions" begin
    m,ps,st=A.load_resolved(joinpath(@__DIR__,"../data/ace_ref/resolved_model_v1.toml"))
    m=A.project_species(m,[:Ti],0:2)
    ch=A.CallableRadial(A.SmoothExponential(),A.RadialMetadata(m.radial_spec;cutoff=5.5,ad=true);
        config=(cutoff=5.5,powers=[n-1 for (n,l) in m.radial_spec]))
    replaced=A.replace_radial(m,ps,fill(ch,2,2);parameters=(rate=[0.4],amplitude=[0.3]))
    m=replaced.model; ps=A.trainable_parameters(m,replaced.parameters)
    sys=TestSystem([0.1 1.6 2.3;0.2 0.4 1.1;0.3 0.7 0.2],Matrix(8.0I,3,3),(false,false,false),[:Ti,:Al,:Ti])
    c=A.ACECalculator(m,ps,st); n=NL.prepare_neighbors(sys,c)
    b=(energy=0.4,forces=reshape(sin.(1:9),3,3),virial=reshape(cos.(1:9),3,3),features=nothing)
    model=A.ScalarSiteModel(m); g=NL.site_pullback(n.positions,model,ps,n,b;centers=[1,3])
    t=PP.parameter_flatten(ps); v=cos.(1:length(t)); eps=2e-6
    obj(t)=contraction(NL.site_evaluate(n.positions,model,PP.parameter_rebuild(ps,t),n;centers=[1,3]),b)
    @test dot(PP.parameter_flatten(g.parameters),v)≈(obj(t+eps*v)-obj(t-eps*v))/(2eps) rtol=1e-5 atol=2e-7
    ref=A.ace_vjp(n.positions,m,ps,n.environment,n.cell,b.energy,b.forces,b.virial;centers=[1,3])
    @test PP.parameter_flatten(g.parameters)≈PP.parameter_flatten(ref.parameters) atol=1e-8
    @test vec(g.positions)≈fd(x -> contraction(NL.site_evaluate(reshape(x,3,3),model,ps,n;centers=[1,3]),b),vec(n.positions)) rtol=2e-5 atol=1e-6
end
struct NoParameterSite end
NI.make_scratch(::NoParameterSite,c,T)=nothing
NI.make_site_output(::NoParameterSite,c,T,o)=NI.make_site_output(NonlinearSite(),c,T,o)
NI.site_derivative_capabilities(::NoParameterSite)=NI.site_derivative_capabilities(NonlinearSite())
NI.evaluate_site!(out,::NoParameterSite,s,ps,scratch,outputs)=NI.evaluate_site!(out,NonlinearSite(),s,(a=0.4,nested=(b=[0.1],)),scratch,outputs)
function NI.site_vjp!(rb,pb,::NoParameterSite,s,ps,bars,scratch;request=NI.SiteInputs())
    NI.site_vjp!(rb,nothing,NonlinearSite(),s,(a=0.4,nested=(b=[0.1],)),bars,scratch;request=NI.SiteInputs(request.displacements,false))
end
@testset "forced splitting and empty parameter derivatives" begin
    sys=TestSystem([3.3 0.0 1.1 2.2;zeros(2,4)],Matrix(8.0I,3,3),(false,false,false),fill(:H,4))
    n=NL.prepare_neighbors(sys,1.3;width=4)
    n=NL.ChunkNeighbours(n.environment,n.positions,n.cell,n.p,n.inverse,[1:4],n.version)
    @test length(NL.form_numa_workers(n,1;scatter_capacity=3).chunks)>1
    p=NoParameterSite(); ps=NamedTuple()
    outputs=NI.SiteOutputs(true,true,true)
    b=(energy=0.2,forces=reshape(sin.(1:12),3,4),virial=reshape(cos.(1:9),3,3),features=nothing)
    g=NL.site_pullback(n.positions,p,ps,n,b;outputs,scatter_capacity=3,batch_size=3)
    @test g.parameters==NamedTuple()
    @test vec(g.positions)≈fd(x -> contraction(NL.site_evaluate(reshape(x,3,4),p,ps,n;outputs,scatter_capacity=3,batch_size=3),b),vec(n.positions)) rtol=1e-5 atol=1e-8
    @test_throws ArgumentError NL.site_evaluate(n.positions,p,ps,n;scatter_capacity=2)
end

@testset "Float32 scalar numerical chunks" begin
    m,ps,st=A.load_resolved(joinpath(@__DIR__,"../data/ace_ref/resolved_model_v1.toml"))
    bundle=A.adapt_model(m,A.trainable_parameters(m,ps);precision=Float32)
    m=bundle.model; ps=bundle.parameters
    sys=TestSystem(Float32[0.1 1.6;0.2 0.4;0.3 0.7],Matrix(8.0f0I,3,3),(false,false,false),[:Ti,:Al])
    c=A.ACECalculator(m,ps,st;precision=Float32); n=NL.prepare_neighbors(sys,c)
    p=A.ScalarSiteModel(m)
    r=NL.site_evaluate(n.positions,p,ps,n)
    ref=A.ace_evaluate(n.positions,m,ps,n.environment,n.cell)
    @test eltype(r.forces)===Float32
    @test r.energy≈ref.energy rtol=1e-6
    @test r.forces≈ref.forces rtol=2e-5
    b=(energy=0.3f0,forces=ones(Float32,3,2),virial=reshape(sin.(Float32.(1:9)),3,3),features=nothing)
    g=NL.site_pullback(n.positions,p,ps,n,b)
    rg=A.ace_vjp(n.positions,m,ps,n.environment,n.cell,b.energy,b.forces,b.virial)
    @test g.positions≈rg.positions rtol=1e-5
    @test PP.parameter_flatten(g.parameters)≈PP.parameter_flatten(rg.parameters) rtol=1e-5
end

@testset "scalar periodic image pullbacks" begin
    m,ps,st=A.load_resolved(joinpath(@__DIR__,"../data/ace_ref/resolved_model_v1.toml"))
    ps=A.trainable_parameters(m,ps)
    sys=TestSystem([0.1 1.2;0.2 0.4;0.3 0.7],Matrix(3.0I,3,3),(true,true,true),[:Ti,:Al])
    c=A.ACECalculator(m,ps,st); n=NL.prepare_neighbors(sys,c)
    b=(energy=0.2,forces=reshape(sin.(1:6),3,2),virial=reshape(cos.(1:9),3,3),features=nothing)
    g=NL.site_pullback(n.positions,A.ScalarSiteModel(m),ps,n,b;centers=[1],scatter_capacity=2)
    ref=A.ace_vjp(n.positions,m,ps,n.environment,n.cell,b.energy,b.forces,b.virial;centers=[1])
    @test g.positions≈ref.positions atol=1e-8
    @test PP.parameter_flatten(g.parameters)≈PP.parameter_flatten(ref.parameters) atol=1e-8
end

end # module SitePullbacksTests
