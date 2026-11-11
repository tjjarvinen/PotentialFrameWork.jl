# Cached live-oracle regression; no ACEpotentials/ETACE runtime dependency.
using SpheriCart, TOML, LinearAlgebra, StaticArrays, ForwardDiff
function may30_agrees(y,ref)
    @test ustrip(u"eV",y.energy) ≈ ref.energy atol=2e-8 rtol=0
    @test reduce(hcat,[ustrip.(u"eV/Å",v) for v in y.forces]) ≈ ref.forces atol=2e-9 rtol=2e-11
    @test ustrip.(u"eV",y.virial) ≈ ref.virial atol=2e-8 rtol=2e-11
end
@testset "May30 solid harmonics and modern pair envelope import" begin
    A=PotentialFrameWork.ACE
    root=joinpath(@__DIR__,"../data/may30")
    path=joinpath(root,"resolved.toml")
    model,ps,st=A.load_resolved(path)
    @test model.angular_kind===:real_solid
    @test model.atomic_numbers==[22,13,23]
    @test size(ps.WB)==(2155,3) && size(ps.Wpair)==(10,3)
    @test first(model.pair_radial).envelope===:PolyEnvelope1sR
    @test length(first(model.radial).nodes)==1000
    @test_throws ArgumentError A.pack_model(model,Array)
    @test_throws ArgumentError A.pack_radials(model.pair_radial,Array,Float64)
    @test_throws ArgumentError A.ACEFeatureLayer(model)
    # Solid = r^l * spherical, including the Cartesian product derivative.
    Rs=[r*SVector(.3,-.7,.9) for r in (.3,1.,2.7)]
    Y,dY=A.angular_ed(A.angular_basis(model,Float64),Rs)
    Ys,dYs=A.angular_ed(A.angular_basis(model.lmax,Float64),Rs)
    for j in eachindex(Rs),l=0:model.lmax,k=l*l+1:(l+1)^2
        r=norm(Rs[j]); g=l==0 ? zero(Rs[j]) : l*r^(l-2)*Rs[j]
        @test Y[j,k] ≈ r^l*Ys[j,k] atol=2e-11 rtol=2e-12
        @test dY[j,k] ≈ r^l*dYs[j,k]+g*Ys[j,k] atol=2e-11 rtol=2e-12
    end
    for T in (Float32,Float64),W in (1,4,8),L in (0,1,6)
        rs=SVector{3,T}.(Rs); b=A.BlockAngular(L;width=W,precision=T,solid=true)
        y=zeros(T,length(rs),(L+1)^2); dy=zeros(SVector{3,T},size(y))
        A.angular_ed!(y,dy,b,rs)
        tol=T===Float64 ? 2e-11 : 4e-4
        @test y ≈ Y[:,1:(L+1)^2] atol=tol rtol=tol
        @test dy ≈ dY[:,1:(L+1)^2] atol=tol rtol=tol
    end
    for ch in model.pair_radial,r in (1.2,2.8,5.99,6.0 - 1e-7,6.,6.0 + 1e-7)
        v,d=A.radial_ed(ch,r)
        @test d ≈ ForwardDiff.derivative(t->A.radial_values(ch,t),r) atol=1e-10 rtol=1e-10
        rc,p=ch.envelope_parameters
        expected=r>=rc ? 0. : ((r/rc)^(-Int(p))-1)*(1-r/rc)
        @test A.envelope(ch,r,A.transformed(ch.transform,r)) ≈ expected
        for W in (1,4,8)
            R=zeros(W,length(v)); dR=similar(R)
            A.radial_block!(R,dR,ntuple(_->r,W),ch,1,Val(W),ntuple(_->true,W))
            @test R[1,:] ≈ v atol=1e-12 rtol=1e-12
            @test dR[1,:] ≈ d atol=1e-12 rtol=1e-12
        end
        if r>=6; @test all(iszero,v) && all(iszero,d); end
    end
    ch=first(model.pair_radial)
    @test_throws DomainError A.radial_values(ch,0.)
    topology(e)=[Set((e.neighbours[k],e.shifts[k]) for k=e.offsets[1,i]:e.offsets[end,i]-1 if e.active[k]) for i in eachindex(e.species)]
    for name in ("bcc-128","hcp-150","omega-240")
        d=TOML.parsefile(joinpath(root,name*".toml")); n=length(d["species"])
        X=reshape(d["positions"],3,n); H=reshape(d["cell"],3,3)
        sys=TestSystem(X,H,(true,true,true),ChemicalSpecies.(d["species"]))
        ref=(;energy=d["energy"],forces=reshape(d["forces"],3,n),virial=reshape(d["virial"],3,3))
        inp=A.ace_inputs(sys,model)
        for algorithm in (:optimized,:dag),backend in (:sphericart,:generated,:block)
            calc=A.ACECalculator(path;slots=1,product_algorithm=algorithm,angular_backend=backend)
            may30_agrees(AC.energy_forces_virial(sys,calc),ref)
        end
        # Functional path and scalar prepared radial path also honor convention.
        functional=A.ace_evaluate(X,model,ps,inp.environment,H)
        @test functional.energy ≈ ref.energy atol=2e-8 rtol=0
        @test functional.forces ≈ ref.forces atol=2e-9 rtol=2e-11
        @test functional.virial ≈ ref.virial atol=2e-8 rtol=2e-11
        ex=A.CPUExecutor(model,inp.environment;slots=1,product_algorithm=:dag)
        A.execute_prepared!(ex,X,H,ps;simd=false)
        @test A.prepared_result(ex).forces ≈ ref.forces atol=2e-9 rtol=2e-11
        if name=="bcc-128"
            small=TestSystem(copy(X[:,1:3]),copy(H),(false,false,false),ChemicalSpecies.(d["species"][1:3]))
            smallinput=A.ace_inputs(small,model)
            sr=A.ace_evaluate(small.xyz,model,ps,smallinput.environment,H)
            gradient=ForwardDiff.gradient(vec(small.xyz)) do x
                A.ace_evaluate(reshape(x,3,:),model,ps,smallinput.environment,H;energy_only=true).energy
            end
            @test reshape(gradient,3,:) ≈ -sr.forces atol=2e-9 rtol=2e-11
            # Verify the solid convention through the existing force-loss VJP.
            B=reshape(sin.(1:9),3,3); direction=reshape(cos.(1:9),3,3)
            vjp=A.ace_vjp(small.xyz,model,ps,smallinput.environment,H,0.,B,zeros(3,3))
            loss(x,p)=sum(B.*A.ace_evaluate(x,model,p,smallinput.environment,H).forces)
            step=2e-5
            numerical=(loss(small.xyz+step*direction,ps)-loss(small.xyz-step*direction,ps))/(2step)
            @test sum(vjp.positions.*direction) ≈ numerical atol=2e-7 rtol=2e-6
            pp=deepcopy(ps); pm=deepcopy(ps)
            pp.WB[10,3]+=step; pm.WB[10,3]-=step
            pp.Wpair[2,3]+=.3step; pm.Wpair[2,3]-=.3step
            numerical=(loss(small.xyz,pp)-loss(small.xyz,pm))/(2step)
            @test vjp.parameters.WB[10,3]+.3vjp.parameters.Wpair[2,3] ≈ numerical atol=2e-7 rtol=2e-6
            energy(x,h)=begin
                A.execute_prepared!(ex,x,h,ps;simd=true)
                A.prepared_result(ex).energy
            end
            h=2e-5; xp=copy(X); xm=copy(X); xp[2,1]+=h; xm[2,1]-=h
            @test -(energy(xp,H)-energy(xm,H))/(2h) ≈ ref.forces[2,1] atol=3e-6 rtol=1e-6
            for x in (xp,xm)
                fresh=A.ace_inputs(TestSystem(x,H,(true,true,true),sys.labels),model).environment
                @test topology(fresh)==topology(inp.environment)
            end
            strain=[.2 .3 -.1;-.17 -.09 .12;.04 -.08 .15]
            plus=I+h*strain; minus=I-h*strain
            @test (energy(plus*X,plus*H)-energy(minus*X,minus*H))/(2h) ≈ -sum(ref.virial.*transpose(strain)) atol=5e-6 rtol=1e-6
            for transform in (plus,minus)
                fresh=A.ace_inputs(TestSystem(transform*X,transform*H,(true,true,true),sys.labels),model).environment
                @test topology(fresh)==topology(inp.environment)
            end
            mktempdir() do dir
                native=joinpath(dir,"fit.ace")
                A.save_calculator(native,A.ACECalculator(path;slots=1))
                restored=A.ACECalculator(native;slots=1)
                @test restored.model.angular_kind===:real_solid
                may30_agrees(AC.energy_forces_virial(sys,restored),ref)
            end
            for T in (Float32,Float64)
                b=A.adapt_model(model,ps;precision=T)
                @test b.model.angular_kind===:real_solid
                @test A.model_with(b.model).angular_kind===:real_solid
            end
        end
    end
end
