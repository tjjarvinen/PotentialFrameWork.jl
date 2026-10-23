using Zygote, ForwardDiff, ChainRulesCore, LuxCore

@testset "prepared ACE SIMD, shared forces and slot schedules" begin
    model,ps,st=ACECore.load_resolved(ACE_MODEL_PATH)
    H=[3.7 0.7 0.1;0 3.8 0.3;0 0 3.6]
    sys=TestSystem([0.1 2.3 1.1;0.2 0.7 2.4;0.3 1.2 1.5],H,(true,true,true),[:Ti,:Al,:Ti])
    for width in (1,2,4,8,16), T in (Float32,Float64)
        inp=ACECore.ace_inputs(sys,model;environment=ACECore.FullEnvironment(;width))
        env=inp.environment
        reference=ACECore.ace_evaluate(inp.positions,model,ps,env,inp.cell)
        w=ACECore.CPUExecutor(model,env;precision=T,slots=3)
        for simd in (false,true), parallel in (false,true)
            ACECore.execute_prepared!(w,inp.positions,inp.cell,ps;simd,parallel)
            r=ACECore.prepared_result(w)
            tol=T===Float64 ? 2e-10 : 3e-4
            @test r.energy ≈ reference.energy atol=tol rtol=tol
            @test r.forces ≈ reference.forces atol=tol rtol=tol
            @test r.virial ≈ reference.virial atol=tol rtol=tol
            @test sizeof(w.force)==3length(sys)*sizeof(T)
        end
        # Every observed site destination belongs to the advertised footprint.
        for i in eachindex(env.species)
            fill!(w.force,0); fill!(w.virials,0)
            ACECore.cpu_site!(w,w.scratch[1],i,inp.positions,inp.cell,ps,Val(width),true)
            actual=[w.p[j] for j=1:length(sys) if any(!iszero,w.force[:,j])]
            @test issubset(actual,env.footprints[i])
        end
        saved=ACECore.prepared_result(w)
        ACECore.prepare_environment!(env,sys,model)
        @test_throws ArgumentError ACECore.execute_prepared!(w,inp.positions,inp.cell,ps)
        ACECore.prepare_executor!(w,env)
        updated=ACECore.install_readout(ps,0.7ACECore.linear_parameters(ps))
        ACECore.execute_prepared!(w,inp.positions,inp.cell,updated)
        expected=ACECore.ace_evaluate(inp.positions,model,updated,env,inp.cell)
        @test ACECore.prepared_result(w).forces ≈ expected.forces atol=3e-4 rtol=3e-4
        @test saved.forces!==w.force
        # Shrink then grow: all active regions reset, including empty groups.
        small=TestSystem(zeros(3,1),H,(false,false,false),[:Ti])
        si=ACECore.ace_inputs(small,model;environment=env)
        ACECore.prepare_executor!(w,env); ACECore.execute_prepared!(w,si.positions,si.cell,ps)
        @test all(iszero,ACECore.prepared_result(w).forces)
        @test ACECore.prepared_result(w).energy ≈ model.reference[1] atol=1e-3
    end
    # A phase with several disjoint sites exercises concurrent writes, rather
    # than relying only on dense tiny fixtures whose phases contain one centre.
    separated=TestSystem([0. 2.5 15. 17.5 30. 32.5;zeros(2,6)],50Matrix{Float64}(I,3,3),
        (false,false,false),[:Ti,:Al,:Ti,:Al,:Ti,:Al])
    inp=ACECore.ace_inputs(separated,model;environment=ACECore.FullEnvironment(width=4))
    @test maximum(length,inp.environment.phases)>1
    w=ACECore.CPUExecutor(model,inp.environment;slots=3)
    serial=ACECore.ace_evaluate(inp.positions,model,ps,inp.environment,inp.cell)
    for _=1:4
        ACECore.execute_prepared!(w,inp.positions,inp.cell,ps;parallel=true)
        @test ACECore.prepared_result(w).forces ≈ serial.forces atol=1e-11
    end
end

@testset "ACE energy/force/virial reverse products and Lux" begin
    model,readout,st=ACECore.load_resolved(ACE_MODEL_PATH)
    ps=ACECore.trainable_parameters(model,readout)
    H=Matrix{Float64}(I,3,3)*9
    sys=TestSystem([0.1 2.3 1.1;0.2 0.7 2.4;0.3 1.2 1.5],H,(false,false,false),[:Ti,:Al,:Ti])
    inp=ACECore.ace_inputs(sys,model); X=inp.positions; env=inp.environment
    flat=PP.parameter_flatten(ps)
    B=sin.(X); C=reshape(cos.(1:9),3,3)
    y=ACECore.ace_evaluate(X,model,ps,env,H)
    @test Zygote.gradient(x -> ACECore.ace_evaluate(x,model,ps,env,H).energy,X)[1] ≈ -y.forces atol=1e-10
    padded=ACECore.ace_inputs(sys,model;environment=ACECore.FullEnvironment(width=16))
    @test Zygote.gradient(x -> sum(B.*ACECore.ace_evaluate(x,model,ps,padded.environment,H).forces),X)[1] ≈
        Zygote.gradient(x -> sum(B.*ACECore.ace_evaluate(x,model,ps,env,H).forces),X)[1] atol=1e-10
    for kind in (:energy,:forces,:virial,:combined)
        function loss(x,parameters)
            r=ACECore.ace_evaluate(x,model,parameters,env,H)
            kind===:energy && return r.energy
            kind===:forces && return sum(B.*r.forces)
            kind===:virial && return sum(C.*r.virial)
            0.7r.energy+sum(B.*r.forces)+sum(C.*r.virial)
        end
        gx,gp=Zygote.gradient(loss,X,ps)
        @test gx ≈ central_gradient(x -> loss(x,ps),X) atol=3e-6 rtol=3e-6
        @test gx ≈ ForwardDiff.gradient(x -> loss(x,ps),X) atol=1e-9 rtol=1e-9
        # Independent parameter direction checks in every parameter subtree.
        for key in keys(ps)
            direction=PP.parameter_rebuild(ps,zeros(length(flat)))
            target=direction[key] isa NamedTuple ? direction[key].scale : direction[key]
            target .= reshape(sin.(1:length(target)),size(target))
            d=PP.parameter_flatten(direction); h=2e-6
            fd=(loss(X,PP.parameter_rebuild(ps,flat+h*d))-loss(X,PP.parameter_rebuild(ps,flat-h*d)))/(2h)
            @test dot(PP.parameter_flatten(gp),d) ≈ fd atol=3e-6 rtol=3e-6
        end
    end
    _,pb=rrule(ACECore.ace_evaluate,X,model,ps,env,H)
    delta=(;energy=0.3,forces=B,virial=C)
    expected=pb(delta)
    X .+= 0.1; fill!(ps.radial.scale,2); fill!(ps.WB,0)
    ACECore.prepare_environment!(env,sys,model)
    again=pb(delta)
    @test again[2] ≈ expected[2]
    @test again[4].WB ≈ expected[4].WB
    @test again[4].radial.scale ≈ expected[4].radial.scale
    @test pb(ZeroTangent())[2] isa ZeroTangent
    @test pb((;energy=Thunk(() -> 1.),forces=ZeroTangent(),virial=ZeroTangent()))[2] isa AbstractMatrix
    layer=lux_layer(model;parameters=readout)
    initial,state=LuxCore.setup(LuxCore.Random.Xoshiro(1234),layer)
    output,newstate=LuxCore.apply(layer,inp,initial,state)
    @test newstate===state
    @test output.energy ≈ ACECore.ace_evaluate(inp.positions,model,initial,inp.environment,inp.cell).energy
    @test all(isfinite,PP.parameter_flatten(Zygote.gradient(p -> sum(B.*first(LuxCore.apply(layer,inp,p,state)).forces),initial)[1]))
    # Nonzero self-images have zero position derivatives, yet nonzero virial
    # parameter derivatives. Both endpoints must scatter to the same atom.
    selfsys=TestSystem(zeros(3,1),3.7Matrix{Float64}(I,3,3),(true,true,true),[:Ti])
    selfinp=ACECore.ace_inputs(selfsys,model)
    @test all(iszero,Zygote.gradient(x -> sum(ACECore.ace_evaluate(x,model,initial,selfinp.environment,selfinp.cell).forces),selfinp.positions)[1])
    gradient=Zygote.gradient(p -> sum(C.*ACECore.ace_evaluate(selfinp.positions,model,p,selfinp.environment,selfinp.cell).virial),initial)[1]
    direction=reshape(sin.(1:length(initial.WB)),size(initial.WB)); h=2e-6
    f=a -> sum(C.*ACECore.ace_evaluate(selfinp.positions,model,merge(initial,(WB=a,)),selfinp.environment,selfinp.cell).virial)
    @test dot(gradient.WB,direction) ≈ (f(initial.WB+h*direction)-f(initial.WB-h*direction))/(2h) atol=3e-6 rtol=3e-6
end

@testset "ACE short deterministic rollout" begin
    model,ps,st=ACECore.load_resolved(ACE_MODEL_PATH)
    H=Matrix{Float64}(I,3,3)*9
    sys=TestSystem([0.1 2.5;0.2 0.7;0.3 1.2],H,(false,false,false),[:Ti,:Al])
    inp=ACECore.ace_inputs(sys,model); X=inp.positions; velocity=0.01cos.(X)
    dt=1e-4
    function rollout(x,v,weights)
        params=merge(ps,(;WB=weights))
        for _=1:2
            f=ACECore.ace_evaluate(x,model,params,inp.environment,H).forces
            vh=v+(dt/2)*f; x=x+dt*vh
            v=vh+(dt/2)*ACECore.ace_evaluate(x,model,params,inp.environment,H).forces
        end
        sum(abs2,x)+0.3sum(abs2,v)
    end
    gx,gv,gp=Zygote.gradient(rollout,X,velocity,ps.WB)
    @test gx ≈ central_gradient(x -> rollout(x,velocity,ps.WB),X) atol=3e-7
    @test gv ≈ central_gradient(v -> rollout(X,v,ps.WB),velocity) atol=3e-7
    direction=reshape(sin.(1:length(ps.WB)),size(ps.WB)); h=2e-6
    @test dot(gp,direction) ≈ (rollout(X,velocity,ps.WB+h*direction)-rollout(X,velocity,ps.WB-h*direction))/(2h) atol=3e-7
end
