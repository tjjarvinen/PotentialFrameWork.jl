using LinearAlgebra, SparseArrays, StaticArrays, SpheriCart, ChainRulesCore, ForwardDiff, Zygote, LuxCore

function feature_central_gradient(f,x;h=2e-6)
    g=similar(x)
    for i in eachindex(x)
        plus=copy(x); minus=copy(x); plus[i]+=h; minus[i]-=h
        g[i]=(f(plus)-f(minus))/(2h)
    end
    g
end
feature_tangent_flatten(x,template::AbstractArray)=vec(x)
feature_tangent_flatten(x,template::NamedTuple)=reduce(vcat,
    (feature_tangent_flatten(getproperty(x,k),getproperty(template,k)) for k in keys(template));init=Float64[])

function feature_radial(;lmax=1)
    spec=lmax==1 ? [(1,0),(1,1),(2,1)] : [(1,0),(1,1),(2,1),(1,2)]
    ACECore.CallableRadial(ACECore.SmoothExponential(),ACECore.RadialMetadata(spec;cutoff=4.,ad=true);
        config=(cutoff=4.,powers=lmax==1 ? [0,0,1] : [0,0,1,0]))
end
feature_radial_ps()=(rate=[0.4],amplitude=[1.3])
function feature_fixture(;correlation=2,lmax=1,species=[:Si,:O],extra=false)
    irreps=((L=0,parity=:even,channels=3),(L=1,parity=:odd,channels=2),(L=1,parity=:even,channels=2))
    extra && (irreps=(irreps...,(L=3,parity=:odd,channels=2),(L=2,parity=:even,channels=2)))
    layer=ACECore.ACEFeatureLayer(;species,radial=feature_radial(;lmax),correlation,lmax,output_irreps=irreps)
    ps=ACECore.feature_parameters(layer;custom_radial=feature_radial_ps())
    # Activate all basis channels, not only the default leading identity rows.
    for (b,W) in enumerate(values(ps.mixing)); W .= reshape(sin.(b .+ (1:length(W))),size(W)); end
    layer,ps
end

@testset "equivariant feature constructor and basis" begin
    layer,ps=feature_fixture()
    summary=ACECore.feature_basis_summary(layer)
    @test [(s.L,s.parity) for s in summary]==[(0,:even),(1,:odd),(1,:even)]
    @test all(s.basis_channels>0 for s in summary)
    @test maximum(length,layer.products)==2
    @test any(d[1]==2 for d in layer.density) # explicit neighbor-species channels
    @test ACECore.feature_basis_summary(first(feature_fixture()))==summary
    @test first(feature_fixture()).coupling==layer.coupling
    base=(species=[:Si],radial=feature_radial(),correlation=2,lmax=1)
    for spec in (((L=-1,parity=:even,channels=1),),((L=1.0,parity=:odd,channels=1),),
                 ((L=true,parity=:odd,channels=1),),((L=0,parity=:even,channels=0),),
                 ((L=0,parity=:bad,channels=1),),((L=0,parity=:even,channels=1),(L=0,parity=:even,channels=2)),
                 ((L=0,parity=:odd,channels=1),),((L=5,parity=:even,channels=1),),())
        @test_throws ArgumentError ACECore.ACEFeatureLayer(;base...,output_irreps=spec)
    end
    @test_throws ArgumentError ACECore.ACEFeatureLayer(;base...,correlation=0,output_irreps=((L=0,parity=:even,channels=1),))
    @test_throws ArgumentError ACECore.ACEFeatureLayer(;base...,cutoff=1.,output_irreps=((L=0,parity=:even,channels=1),))
    @test_throws ArgumentError ACECore.ACEFeatureLayer(;base...,totaldegree=1,output_irreps=((L=1,parity=:odd,channels=1),))
    @test_throws ArgumentError ACECore.ACEFeatureLayer(;base...,species=[:Si,:Si],output_irreps=((L=0,parity=:even,channels=1),))
    # A single 1o density cannot produce an axial 1e by self coupling.
    ch=ACECore.CallableRadial(ACECore.SmoothExponential(),ACECore.RadialMetadata([(1,1)];cutoff=4.,ad=true);config=(cutoff=4.,powers=[0]))
    @test_throws ArgumentError ACECore.ACEFeatureLayer(species=[:Si],radial=ch,correlation=2,lmax=1,output_irreps=((L=1,parity=:even,channels=1),))
    symmetric=ACECore.ACEFeatureLayer(species=[:Si],radial=ch,correlation=3,lmax=1,
        output_irreps=((L=1,parity=:odd,channels=1),(L=3,parity=:odd,channels=1),(L=2,parity=:even,channels=1)))
    @test [s.basis_channels for s in ACECore.feature_basis_summary(symmetric)]==[2,1,1]
    wide=ACECore.ACEFeatureLayer(;base...,output_irreps=((L=1,parity=:even,channels=8),))
    @test only(ACECore.feature_basis_summary(wide)).basis_channels==1
    @test only(ACECore.feature_basis_summary(wide)).rank_limit==1
end

@testset "resolved cross product reference" begin
    # Small explicit product table in the SpheriCart (y,z,x) convention:
    # 0e=A_0, 1o=A_1, 1e=-cross(A_1,A_2)/sqrt(2).
    radial=fill(feature_radial(),1,1)
    density=[(1,1,1);[(1,n,m) for n in 2:3 for m in 2:4]]
    products=[[a] for a=1:4]
    crossproducts=[[2,6],[3,7],[4,5],[2,7],[3,5],[4,6]]
    append!(products,crossproducts)
    C0=sparse([1],[1],[1.],1,10)
    C1=sparse(1:3,2:4,ones(3),3,10)
    # -cross in y,z,x order: (x1*z2-z1*x2, y1*x2-x1*y2, z1*y2-y1*z2)
    C2=sparse([1,1,2,2,3,3],[10,6,8,7,9,5],[1.,-1,1,-1,1,-1]/sqrt(2),3,10)
    specs=ACECore.feature_irreps(((L=0,parity=:even,channels=1),(L=1,parity=:odd,channels=1),(L=1,parity=:even,channels=1)))
    l=ACECore.ACEFeatureLayer(radial,density,products,(C0,C1,C2);species=[:Si],radial_spec=[(1,0),(1,1),(2,1)],lmax=1,output_irreps=((L=0,parity=:even,channels=1),(L=1,parity=:odd,channels=1),(L=1,parity=:even,channels=1)))
    X=[0.1 1.2 0.7 2.1;0.2 0.4 1.8 1.1;0.3 1.1 0.5 1.9]; H=9Matrix{Float64}(I,3,3)
    inp=ACECore.ace_inputs(TestSystem(X,H,(false,false,false),fill(:Si,4)),l)
    ps=ACECore.feature_parameters(l;custom_radial=feature_radial_ps())
    y=ACECore.ace_features(X,l,ps,inp.environment,inp.cell)
    generated=ACECore.ACEFeatureLayer(species=[:Si],radial=feature_radial(),correlation=2,lmax=1,output_irreps=((L=1,parity=:even,channels=1),))
    pg=ACECore.feature_parameters(generated;custom_radial=feature_radial_ps())
    yg=only(ACECore.ace_features(X,generated,pg,inp.environment,inp.cell))
    for i=1:4
        a=zeros(3); b=zeros(3); scalar=0.
        for j=1:4
            i==j && continue
            d=X[:,j]-X[:,i]; r=norm(d)
            v=ACECore.radial_values(radial[1,1],r,ps.custom_radial)
            scalar+=v[1]/sqrt(4pi)
            a+=v[2]*sqrt(3/4pi)*d/r; b+=v[3]*sqrt(3/4pi)*d/r
        end
        expected=-cross(a,b)[[2,3,1]]/sqrt(2)
        @test y[1][1,1,i] ≈ scalar
        @test y[2][:,1,i] ≈ a[[2,3,1]]
        @test y[3][:,1,i] ≈ expected atol=1e-15
        @test yg[:,1,i] ≈ expected atol=1e-15
        @test norm(expected)>1e-5
    end
end

# Obtain D_L(R) independently from sampled SpheriCart harmonics, rather than
# reusing the coupling machinery under test. D acts on column m components.
function feature_rotation(L,R)
    basis=ACECore.angular_basis(L,Float64)
    points=[normalize(SVector(sin(0.7i),cos(1.3i),sin(2.1i+0.2))) for i=1:2L+7]
    Y=first(ACECore.angular_ed(basis,points))[:,L^2+1:(L+1)^2]
    Yrot=first(ACECore.angular_ed(basis,[SVector{3}(R*p) for p in points]))[:,L^2+1:(L+1)^2]
    transpose(Y\Yrot)
end

@testset "generated feature symmetries, correlation three and higher L" begin
    layer,ps=feature_fixture(correlation=3,lmax=2,extra=true)
    @test maximum(length,layer.products)==3
    X=[0.1 1.2 0.7 2.1;0.2 0.4 1.8 1.1;0.3 1.1 0.5 1.9]; H=9Matrix{Float64}(I,3,3)
    labels=[:Si,:O,:Si,:O]; sys=TestSystem(X,H,(false,false,false),labels)
    inp=ACECore.ace_inputs(sys,layer;environment=ACECore.FullEnvironment(width=8)); env=inp.environment
    y=ACECore.ace_features(X,layer,ps,env,inp.cell)
    @test y isa Tuple
    for (s,block) in zip(layer.irreps,y)
        @test size(block)==(2s.L+1,s.channels,4)
        @test norm(block)>1e-5
    end
    axis=normalize([1.,2.,-0.7]); K=[0. -axis[3] axis[2];axis[3] 0. -axis[1];-axis[2] axis[1] 0.]
    R=Matrix{Float64}(I,3,3)+sin(0.73)*K+(1-cos(0.73))*K*K
    rotated=ACECore.ace_features(R*X,layer,ps,env,R*inp.cell)
    inverted=ACECore.ace_features(-X,layer,ps,env,-inp.cell)
    translated=ACECore.ace_features(X .+ [0.3,-0.7,1.4],layer,ps,env,inp.cell)
    perm=[3,1,4,2]; pi=ACECore.ace_inputs(TestSystem(X[:,perm],H,sys.pbc,labels[perm]),layer)
    permuted=ACECore.ace_features(pi.positions,layer,ps,pi.environment,pi.cell)
    for (b,s) in enumerate(layer.irreps)
        D=feature_rotation(s.L,R)
        for i=1:4; @test rotated[b][:,:,i] ≈ D*y[b][:,:,i] atol=2e-12 rtol=2e-10; end
        @test inverted[b] ≈ (s.parity===:even ? 1 : -1)*y[b] atol=1e-13
        @test translated[b] ≈ y[b] atol=1e-13
        @test permuted[b] ≈ y[b][:,:,perm] atol=1e-13
    end
    again=ACECore.ace_features(X,layer,ps,env,inp.cell)
    @test all(again[b]!==y[b] && again[b]==y[b] for b in eachindex(y))
end

@testset "feature pullback positions, parameters, periodic topology and lifetime" begin
    l,ps=feature_fixture()
    X=[0.1 1.2 0.7;0.2 0.4 1.8;0.3 1.1 0.5]; H=[3.8 0.2 0.;0. 4.1 0.1;0. 0. 4.3]
    for pbc in ((false,false,false),(true,true,true))
        inp=ACECore.ace_inputs(TestSystem(X,H,pbc,[:Si,:O,:Si]),l;environment=ACECore.FullEnvironment(width=4))
        env=inp.environment; h=2e-6
        y,pb=rrule(ACECore.ace_features,X,l,ps,env,inp.cell)
        cot=map(b -> reshape(sin.(1:length(b)),size(b)),y)
        for active in (nothing,2,3)
            c=map(eachindex(y)) do b; isnothing(active) || active==b ? cot[b] : zero(cot[b]); end
            c=Tuple(c)
            loss(x,p)=sum(sum(a.*b) for (a,b) in zip(ACECore.ace_features(x,l,p,env,inp.cell),c))
            gx,gp=Zygote.gradient(loss,X,ps)
            @test gx ≈ feature_central_gradient(x -> loss(x,ps),X) atol=2e-8 rtol=3e-6
            @test gx ≈ ForwardDiff.gradient(x -> loss(x,ps),X) atol=1e-11 rtol=1e-10
            @test sum(gx;dims=2) ≈ zeros(3,1) atol=1e-12
            flat=PP.parameter_flatten(ps)
            for key in keys(ps)
                direction=PP.parameter_rebuild(ps,zeros(length(flat)))
                function activate!(tree)
                    tree isa NamedTuple ? foreach(activate!,values(tree)) : (tree .= reshape(cos.(1:length(tree)),size(tree)))
                end
                activate!(getproperty(direction,key)); d=PP.parameter_flatten(direction)
                fd=(loss(X,PP.parameter_rebuild(ps,flat+h*d))-loss(X,PP.parameter_rebuild(ps,flat-h*d)))/(2h)
                @test dot(PP.parameter_flatten(gp),d) ≈ fd atol=2e-8 rtol=3e-6
            end
            thunked=Tuple(isnothing(active) || active==b ? Thunk(() -> cot[b]) : ZeroTangent() for b in eachindex(y))
            @test pb(thunked)[2] ≈ gx atol=1e-12
            @test pb(c)[4].mixing.block_2 ≈ gp.mixing.block_2 atol=1e-12
        end
        expected=pb(cot); saved=deepcopy((X,ps))
        # Caller mutation and a later forward/topology reuse cannot alter saves.
        X .+= 0.03; ps.mixing.block_1 .*= 2; ps.radial.scale .*= 0.8; ps.custom_radial.rate .+= 0.1
        ACECore.ace_features(X,l,ps,env,inp.cell)
        ACECore.prepare_environment!(env,TestSystem(X,H,pbc,[:Si,:O,:Si]),l)
        @test pb(cot)[2] ≈ expected[2]
        @test feature_tangent_flatten(pb(cot)[4],ps) ≈ feature_tangent_flatten(expected[4],ps)
        @test pb(ZeroTangent())[2] isa ZeroTangent
        @test all(iszero,pb(map(zero,y))[2])
        X,ps=saved
    end
    for n in (0,1)
        inp=ACECore.ace_inputs(TestSystem(zeros(3,n),H,(false,false,false),fill(:Si,n)),l)
        y,pb=rrule(ACECore.ace_features,inp.positions,l,ps,inp.environment,inp.cell)
        @test all(all(iszero,b) for b in y)
        g=pb(map(b -> ones(size(b)),y))
        @test all(iszero,g[2])
        @test all(iszero,feature_tangent_flatten(g[4],ps))
    end
    self=ACECore.ace_inputs(TestSystem(zeros(3,1),3.2Matrix{Float64}(I,3,3),(true,true,true),[:Si]),l)
    sy,spb=rrule(ACECore.ace_features,self.positions,l,ps,self.environment,self.cell)
    sg=spb(map(b -> ones(size(b)),sy))
    @test norm(first(sy))>0
    @test all(iszero,sg[2])
    @test norm(sg[4].mixing.block_1)>0
    # A separate invariant head gives conservative forces through this VJP.
    inp=ACECore.ace_inputs(TestSystem(X,H,(false,false,false),[:Si,:O,:Si]),l)
    energy(x,p)=sum(sum(abs2,b) for b in ACECore.ace_features(x,l,p,inp.environment,inp.cell))
    force=-Zygote.gradient(x -> energy(x,ps),X)[1]
    @test force ≈ -feature_central_gradient(x -> energy(x,ps),X) atol=2e-8 rtol=3e-6
    layer=lux_layer(l;parameters=ps)
    initial,state=LuxCore.setup(LuxCore.Random.Xoshiro(7),layer)
    output,newstate=LuxCore.apply(layer,inp,initial,state)
    @test newstate===state
    @test output==ACECore.ace_features(X,l,ps,inp.environment,inp.cell)
    @test all(isfinite,PP.parameter_flatten(Zygote.gradient(p -> sum(first(first(LuxCore.apply(layer,inp,p,state)))) ,initial)[1]))
end

@testset "resolved scalar many-body feature equality" begin
    m,readout,_=ACECore.load_resolved(ACE_MODEL_PATH)
    l=ACECore.ACEFeatureLayer(m); p=ACECore.feature_parameters(l)
    p.mixing.block_1[:,1,:].=readout.WB
    X=[0.1 2.3 1.1;0.2 0.7 2.4;0.3 1.2 1.5]; H=7Matrix{Float64}(I,3,3)
    inp=ACECore.ace_inputs(TestSystem(X,H,(true,true,true),[:Ti,:Al,:Ti]),m)
    E(x)=sum(only(ACECore.ace_features(x,l,p,inp.environment,inp.cell)))
    scalar=ACECore.ace_evaluate(X,m,readout,inp.environment,inp.cell)
    baseline=ACECore.ace_evaluate(X,m,(;WB=zero(readout.WB),Wpair=readout.Wpair),inp.environment,inp.cell)
    @test E(X) ≈ scalar.energy-baseline.energy atol=1e-10
    @test -Zygote.gradient(E,X)[1] ≈ scalar.forces-baseline.forces atol=1e-10
    @test sum(m.reference[inp.environment.species]) != baseline.energy
end
