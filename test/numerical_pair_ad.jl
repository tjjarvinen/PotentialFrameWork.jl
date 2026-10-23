using ForwardDiff, ChainRulesCore, Zygote, LuxCore

function central_gradient(f,x;h=2e-6)
    map(eachindex(x)) do k
        xp=copy(x); xm=copy(x); xp[k]+=h; xm[k]-=h
        (f(xp)-f(xm))/(2h)
    end |> g -> reshape(g,size(x))
end

@testset "pair LuxCore property interface" begin
    # Load the actual optional interface, including structured parameter rules.
    p=potential(); s=fixture(3;cluster=true,pbc=(false,false,false)); inp=pair_inputs(s,p)
    layer=lux_layer(p)
    ps,st=LuxCore.setup(LuxCore.Random.Xoshiro(1234),layer)
    y,returned=LuxCore.apply(layer,inp,ps,st)
    @test returned===st
    @test y.energy ≈ pair_evaluate(inp.positions,p,ps,inp.environment,inp.cell).energy
    grad=Zygote.gradient(a -> sum(inp.positions.*first(LuxCore.apply(layer,inp,(;amplitudes=a),st)).forces),ps.amplitudes)[1]
    @test grad ≈ central_gradient(a -> sum(inp.positions.*pair_evaluate(inp.positions,p,(;amplitudes=a),inp.environment,inp.cell).forces),ps.amplitudes) atol=1e-8
end

@testset "explicit pair geometry and derivative products" begin
    H=[1.7 0.6 0.2;0 1.8 0.4;0 0 1.6]
    small=TestSystem(H*[0.1 0.4 0.8;0.2 0.6 0.3;0.1 0.7 0.5],H,(true,true,true),[:H,:He,:H])
    lj,_=ArgonComparison.argon_models()
    for (sys,p) in ((small,potential()),(fixture(5;cluster=true),potential()),
                    (ArgonComparison.argon_pair(4.2),lj),
                    (fixture(0),potential()),
                    (TestSystem(zeros(3,1),H,(true,true,true),[:H]),potential()))
        inp=pair_inputs(sys,p); X=inp.positions; env=inp.environment; cell=inp.cell
        ps=AC.get_parameters(p)
        r=pair_evaluate(X,p,ps,env,cell)
        agrees((;energy=r.energy*u"eV",forces=r.forces*u"eV/Å",virial=r.virial*u"eV"),
            (;energy=ustrip(evaluate(sys,p).energy),forces=ustrip.(evaluate(sys,p).forces),virial=ustrip.(evaluate(sys,p).virial)))
        B=reshape(sin.(1:length(X)),size(X)); C=reshape(cos.(1:9),3,3)
        energy=x -> pair_evaluate(x,p,ps,env,cell).energy
        @test Zygote.gradient(energy,X)[1] ≈ -r.forces atol=1e-11
        @test ForwardDiff.gradient(energy,X) ≈ -r.forces atol=1e-11
        for kind in (:energy,:forces,:virial,:combined)
            function loss(x,params)
                y=pair_evaluate(x,p,params,env,cell)
                kind===:energy && return y.energy
                kind===:forces && return sum(B.*y.forces)
                kind===:virial && return sum(C.*y.virial)
                0.7y.energy+sum(B.*y.forces)+sum(C.*y.virial)
            end
            gx,gp=Zygote.gradient(loss,X,ps)
            @test gx ≈ central_gradient(x -> loss(x,ps),X) atol=2e-7 rtol=2e-6
            @test gx ≈ ForwardDiff.gradient(x -> loss(x,ps),X) atol=1e-10
            for name in keys(ps)
                a=ps[name]
                if a isa Real
                    f=v -> loss(X,merge(ps,NamedTuple{(name,)}((v[1],))))
                    @test gp[name] ≈ central_gradient(f,[a])[1] atol=2e-7 rtol=2e-6
                else
                    f=v -> loss(X,merge(ps,NamedTuple{(name,)}((v,))))
                    @test gp[name] ≈ central_gradient(f,a) atol=2e-7 rtol=2e-6
                end
            end
        end
        # Saved pullback owns inputs, maps and parameters; repeat after mutation.
        _,pb=rrule(pair_evaluate,X,p,ps,env,cell)
        delta=(;energy=0.7,forces=B,virial=C)
        before=pb(delta)
        X .+= 0.3; cell .*= 1.1
        for a in values(ps); a isa AbstractArray && fill!(a,3); end
        fill!(env.centres,1)
        after=pb(delta)
        @test before[2] ≈ after[2]
        @test before[4] == after[4]
        @test pb(ZeroTangent())[2] isa ZeroTangent
        @test pb((;energy=Thunk(() -> 1.),forces=ZeroTangent(),virial=ZeroTangent()))[2] isa AbstractMatrix
    end
end

@testset "short fixed-cell pair rollout sensitivities" begin
    s=fixture(3;cluster=true,pbc=(false,false,false)); p=potential()
    inp=pair_inputs(s,p); X=inp.positions; env=inp.environment; H=inp.cell
    ps=AC.get_parameters(p); velocity=0.01cos.(X); dt=0.005
    function rollout(x,v,a)
        ps=(;amplitudes=a)
        for _=1:3
            f=pair_evaluate(x,p,ps,env,H).forces
            vh=v+(dt/2)*f
            x=x+dt*vh
            v=vh+(dt/2)*pair_evaluate(x,p,ps,env,H).forces
        end
        sum(abs2,x)+0.3sum(abs2,v)
    end
    gx,gv,gp=Zygote.gradient(rollout,X,velocity,ps.amplitudes)
    @test gx ≈ central_gradient(x -> rollout(x,velocity,ps.amplitudes),X) atol=1e-8
    @test gv ≈ central_gradient(v -> rollout(X,v,ps.amplitudes),velocity) atol=1e-8
    @test gp ≈ central_gradient(a -> rollout(X,velocity,a),ps.amplitudes) atol=1e-8
end
