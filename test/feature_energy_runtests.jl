module FeatureEnergyTests
using PotentialFrameWork, AtomsBase, Test, LinearAlgebra, StaticArrays, SpheriCart, ChainRulesCore, ForwardDiff, Zygote, LuxCore
include("../examples/site_training_support.jl")
const A=PotentialFrameWork.ACE
const NL=PotentialFrameWork.NeighbourList
const NI=PotentialFrameWork.NeighbourInterface
const PP=PotentialFrameWork.PairPotentials
using .SiteTrainingSupport: SmoothInvariantHead
function direction_fd(f,x,v;h=2e-6)
    (f(x+h*v)-f(x-h*v))/(2h)
end
@testset "nonlinear feature energy products" begin
    l,fp,n=SiteTrainingSupport.setup()
    ps=(features=fp,head=(weights=[0.7,-0.3,0.4,0.6],bias=0.2,amplitude=0.8,quadratic=0.11))
    model=A.ACEFeatureEnergy(l,SmoothInvariantHead()); X=n.positions
    y=A.feature_energy_evaluate(X,model,ps,n;batch_size=3)
    E(x,p)=sum(model.head(Tuple(h[:,:,i] for h in A.ace_features(x,l,p.features,n.environment,n.cell)),p.head) for i in axes(x,2))
    @test y.energy≈E(X,ps)
    v=reshape(sin.(1:length(X)),size(X))
    @test -sum(y.forces.*v)≈direction_fd(x -> E(x,ps),X,v) rtol=2e-6 atol=1e-9
    bars=(energy=0.3,forces=reshape(cos.(1:length(X)),size(X)),virial=reshape(sin.(1:9),3,3),features=nothing)
    loss(x,p)=begin
        r=A.feature_energy_evaluate(x,model,p,n)
        bars.energy*r.energy+sum(bars.forces.*r.forces)+sum(bars.virial.*r.virial)
    end
    g=NL.site_pullback(X,A.FeatureEnergySiteModel(model),ps,n,bars)
    flat=PP.parameter_flatten(ps); gp=PP.parameter_flatten(g.parameters)
    # Check each parameter subtree independently, including every mixing block.
    offset=0
    leaves=(values(ps.features.mixing)...,ps.features.radial.scale,
            ps.features.custom_radial.rate,ps.features.custom_radial.amplitude,
            ps.head.weights,ps.head.bias,ps.head.amplitude,ps.head.quadratic)
    for leaf in leaves
        nleaf=leaf isa Number ? 1 : length(leaf)
        d=zeros(length(flat)); d[offset+1:offset+nleaf].=sin.(1:nleaf); offset+=nleaf
        for h in (1e-5,2e-6)
            @test dot(gp,d)≈direction_fd(t -> loss(X,PP.parameter_rebuild(ps,t)),flat,d;h) rtol=1e-5 atol=1e-9
        end
    end
    @test offset==length(flat)
    @test sum(g.positions.*v)≈direction_fd(x -> loss(x,ps),X,v) rtol=1e-5 atol=1e-9
    zg=Zygote.gradient((x,p) -> loss(x,p),X,ps)
    @test zg[1]≈g.positions atol=1e-10
    @test PP.parameter_flatten(zg[2])≈gp atol=1e-10
    h=A.ace_features(X,l,fp,n.environment,n.cell)
    @test all(norm(block)>1e-14 for block in h)
    hb=map(block -> reshape(cos.(1:length(block)),size(block)),h)
    jv=A.feature_jvp(X,l,fp,n.environment,n.cell,v)
    vg=A.feature_vjp(X,l,fp,n.environment,n.cell,hb)
    flatp=PP.parameter_flatten(fp)
    pd=PP.parameter_rebuild(fp,sin.(1:length(flatp)))
    jp=A.feature_jvp(X,l,fp,n.environment,n.cell,zero(X);parameter_direction=pd)
    @test sum(sum(a.*b) for (a,b) in zip(hb,jp))≈dot(PP.parameter_flatten(vg.parameters),PP.parameter_flatten(pd)) atol=1e-10
    @test_throws ArgumentError A.feature_jvp(X,l,fp,n.environment,n.cell,v;parameter_direction=(wrong=ones(length(flatp)),))
    @test sum(sum(a.*b) for (a,b) in zip(hb,jv))≈sum(vg.positions.*v) atol=1e-10
    for b in eachindex(h)
        @test jv[b]≈direction_fd(x -> A.ace_features(x,l,fp,n.environment,n.cell)[b],X,v) rtol=1e-6 atol=1e-9
    end
end

struct NestedRadial end
function (::NestedRadial)(values,derivatives,r,ps,c)
    A.SmoothExponential()(values,derivatives,r,(rate=ps.shape.rate,amplitude=[ps.gain]),c)
end
struct NestedValuesRadial end
(::NestedValuesRadial)(values,r,ps,c)=NestedRadial()(values,nothing,r,ps,c)
function composed_fixture(;periodic=false,values_only=false)
    l,fp,n=SiteTrainingSupport.setup(;periodic)
    old=l.radial[1]
    ch=values_only ? A.CallableRadial(NestedValuesRadial(),old.metadata;config=old.config,derivative=:ad) :
                    A.CallableRadial(NestedRadial(),old.metadata;config=old.config)
    l=A.ACEFeatureLayer(fill(ch,2,2),l.density,l.products,l.coupling;species=[:Si,:O],radial_spec=l.radial_spec,lmax=l.lmax,output_irreps=map(s -> (;L=s.L,parity=s.parity,channels=s.channels),l.irreps))
    fp=merge(fp,(custom_radial=(shape=(rate=[0.4],),gain=1.3),))
    ps=(features=fp,head=(weights=[0.7,-0.3,0.4,0.6],bias=0.2,amplitude=0.8,quadratic=0.11))
    A.ACEFeatureEnergy(l,SmoothInvariantHead()),ps,n
end
@testset "composed periodic lifetime mixed products and nested radials" begin
    for periodic in (false,true), values_only in (false,true)
        model,ps,n=composed_fixture(;periodic,values_only); X=n.positions
        b=(energy=0.1,forces=reshape(cos.(1:9),3,3),virial=reshape(sin.(1:9),3,3),features=nothing)
        opts=(centers=[3,1],batch_size=3,scatter_capacity=3)
        function loss(x,p)
            r=A.feature_energy_evaluate(x,model,p,n;opts...)
            b.energy*r.energy+sum(b.forces.*r.forces)+sum(b.virial.*r.virial)
        end
        g=NL.site_pullback(X,A.FeatureEnergySiteModel(model),ps,n,b;opts...)
        t=PP.parameter_flatten(ps); d=sin.(1:length(t)); v=reshape(cos.(1:9),3,3)
        @test dot(PP.parameter_flatten(g.parameters),d)≈direction_fd(q -> loss(X,PP.parameter_rebuild(ps,q)),t,d) rtol=2e-5 atol=2e-9
        @test sum(g.positions.*v)≈direction_fd(x -> loss(x,ps),X,v) rtol=2e-5 atol=2e-9
        onlyp=NL.site_pullback(X,A.FeatureEnergySiteModel(model),ps,n,b;opts...,request=NI.SiteInputs(false,true))
        @test PP.parameter_flatten(onlyp.parameters)≈PP.parameter_flatten(g.parameters)
        # Arbitrary direct feature bars combine additively with E/F/W bars.
        r=A.feature_energy_evaluate(X,model,ps,n;opts...,outputs=NI.SiteOutputs(true,true,true))
        hb=[map(x -> cos.(x),h) for h in r.features]
        allbars=merge(b,(features=hb,))
        both=NL.site_pullback(X,A.FeatureEnergySiteModel(model),ps,n,allbars;opts...,outputs=NI.SiteOutputs(true,true,true))
        direct=NL.site_pullback(X,A.FeatureSiteModel(model.layer),ps.features,n,(energy=nothing,forces=nothing,virial=nothing,features=hb);opts...,outputs=NI.SiteOutputs(false,false,true))
        @test both.positions≈g.positions+direct.positions
        @test PP.parameter_flatten(both.parameters.features)≈PP.parameter_flatten(g.parameters.features)+PP.parameter_flatten(direct.parameters)
        y,pb=rrule(A.feature_energy_evaluate,X,model,ps,n;opts...)
        saved=pb(b)
        X .+= 0.03; ps.features.custom_radial.shape.rate .+= 0.02
        n.environment.shifts .+= Ref(SVector(1,0,0))
        @test pb(b)[2]≈saved[2]
        @test pb(ZeroTangent())[2] isa AbstractZero
    end
end
@testset "composed symmetry and force-position rollout" begin
    model,ps,n=composed_fixture(); X=n.positions
    p=A.FeatureEnergySiteModel(model)
    B=reshape(sin.(1:9),3,3); C=reshape(cos.(1:9),3,3)
    bars(B)=(energy=nothing,forces=B,virial=nothing,features=nothing)
    gb=NL.site_pullback(X,p,ps,n,bars(B)); gc=NL.site_pullback(X,p,ps,n,bars(C))
    @test sum(C.*gb.positions)≈sum(B.*gc.positions) atol=1e-10
    # Rotate/invert both geometry and force cotangents. Open fixed cell is zero.
    Q=[cos(0.4) -sin(0.4) 0.;sin(0.4) cos(0.4) 0.;0. 0. 1.]
    y=A.feature_energy_evaluate(X,model,ps,n)
    for R in (Q,-Q)
        yr=A.feature_energy_evaluate(R*X,model,ps,n)
        @test yr.energy≈y.energy atol=1e-11
        @test yr.forces≈R*y.forces atol=1e-10
        gr=NL.site_pullback(R*X,p,ps,n,bars(R*B))
        @test gr.positions≈R*gb.positions atol=1e-10
        @test PP.parameter_flatten(gr.parameters)≈PP.parameter_flatten(gb.parameters) atol=1e-10
    end
    function rollout(x,p)
        first=A.feature_energy_evaluate(x,model,p,n)
        next=x+0.01first.forces
        second=A.feature_energy_evaluate(next,model,p,n)
        sum(abs2,second.forces)+0.1second.energy
    end
    gx,gp=Zygote.gradient(rollout,X,ps)
    @test sum(gx.*B)≈direction_fd(x -> rollout(x,ps),X,B) rtol=2e-5 atol=1e-9
    t=PP.parameter_flatten(ps); d=sin.(1:length(t))
    @test dot(PP.parameter_flatten(gp),d)≈direction_fd(q -> rollout(X,PP.parameter_rebuild(ps,q)),t,d) rtol=2e-5 atol=1e-9
end
@testset "composed zero density empty and isolated" begin
    model,ps,n=composed_fixture()
    ps.features.radial.scale .= 0.
    y=A.feature_energy_evaluate(n.positions,model,ps,n)
    @test iszero(y.forces)
    @test iszero(y.virial)
    b=(energy=1.,forces=ones(3,3),virial=ones(3,3),features=nothing)
    g=NL.site_pullback(n.positions,A.FeatureEnergySiteModel(model),ps,n,b)
    @test all(isfinite,PP.parameter_flatten(g.parameters))
    @test all(isfinite,g.positions)
    for count in (0,1)
        sys=SiteTrainingSupport.TrainingSystem(zeros(3,count),Matrix(8.0I,3,3),(false,false,false),fill(AtomsBase.ChemicalSpecies(:Si),count))
        empty_n=NL.prepare_neighbors(sys,3.2;channels=ones(Int,count),species_count=2,width=4)
        r=A.feature_energy_evaluate(empty_n.positions,model,ps,empty_n;batch_size=3)
        @test size(r.forces)==(3,count)
        @test iszero(r.forces)
        gg=Zygote.gradient(p -> A.feature_energy_evaluate(empty_n.positions,model,p,empty_n).energy,ps)[1]
        @test all(isfinite,PP.parameter_flatten(gg))
    end
end

end # module FeatureEnergyTests
