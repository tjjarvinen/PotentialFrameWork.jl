using ForwardDiff, ChainRulesCore, Zygote

struct ChangingCutoff <: AbstractPairPotential
    radius::Float64
end
PotentialFrameWork.cutoff(p::ChangingCutoff)=p.radius
PotentialFrameWork.with_parameters(::ChangingCutoff,ps)=ChangingCutoff(ps.radius)


# This file follows lennard_jones.jl, which loads the shared Argon model.
@testset "Explicit energy parameters and ForwardDiff" begin
    s=fixture(43;L=4.1,pbc=(true,false,true),cluster=true)
    p=potential(); params=AC.get_parameters(p); theta=params.amplitudes
    @test theta==[1.,0.7,1.3]
    @test with_parameters(p,params)!==p
    @test AC.set_parameters!(p,params).amplitudes==p.amplitudes
    @test p.amplitudes==[1. 0.7;0.7 1.3]
    expected=ustrip(u"eV",AC.potential_energy(s,p))
    oracle_gradient=[oracle(s,with_parameters(p,(;amplitudes=Float64.(1:3 .== k)))).energy for k=1:3]
    for width in (1,4,8), parallel in (false,true)
        w=Workspace(;width)
        f=x -> ustrip(u"eV",AC.calculate(AC.Energy(),s,p;
            parameters=(;amplitudes=x),neighbourlist=w,parallel).energy)
        @test f(theta) ≈ expected
        @test ForwardDiff.gradient(f,theta) ≈ oracle_gradient
        @test ForwardDiff.hessian(f,theta) ≈ zeros(3,3) atol=1e-13
        @test validate_schedule(w)
        @test size(w.force)==(3,length(s))
        @test eltype(w.force)===Float64
        @test energy_parameter_gradient(s,p,params;neighbourlist=w,parallel).amplitudes ≈ oracle_gradient
    end
    state=(token=3,)
    for r in (AC.calculate(AC.Energy(),s,p,params,state),
              AC.calculate((AC.Energy(),),s,p,params;state),
              AC.calculate(AC.Energy(),s,p;parameters=params,state))
        @test ustrip(r.energy) ≈ expected
        @test r.state===state
    end
    @test_throws ArgumentError AC.potential_energy(s,p;parameters=params)
    @test_throws ArgumentError AC.forces(s,p;parameters=params)
    @test_throws ArgumentError AC.calculate(AC.Forces(),s,p;parameters=params)
    @test_throws ArgumentError AC.calculate((AC.Energy(),AC.Virial()),s,p;parameters=params)
    @test_throws ArgumentError AC.calculate(AC.Energy(),s,p;parameters=theta)
    @test_throws ArgumentError AC.calculate(AC.Energy(),s,p;parameters=(amplitudes=[1.],))
    @test_throws ArgumentError AC.calculate(AC.Energy(),s,p;parameters=(amplitudes=[1.,NaN,2.],))
    @test_throws ArgumentError AC.calculate(AC.Energy(),s,p;parameters=(radius=3.,))
    @test_throws ArgumentError AC.calculate(AC.Energy(),s,EnergyOnlyPair();parameters=(a=1.,))
    @test_throws ArgumentError parameter_energy(s,ChangingCutoff(2.),(radius=2.1,))
    for n in (0,1)
        empty_or_single=fixture(n)
        f=x -> ustrip(AC.calculate(AC.Energy(),empty_or_single,p;(parameters=(;amplitudes=x),)...).energy)
        @test f(theta)==0
        @test ForwardDiff.gradient(f,theta)==zeros(3)
    end
    # A parameter type change must not round the fixed physical cutoff.
    p2=SmoothPair([:H,:He],p.amplitudes,2.0123456789)
    theta32=Float32.(theta)
    @test cutoff(with_parameters(p2,(;amplitudes=theta32)))===cutoff(p2)
    @test all(isfinite,ForwardDiff.gradient(x -> parameter_energy(s,p2,(;amplitudes=x)),theta32))
    # Geometry precision is independent of parameter precision.
    w=Workspace(width=8,precision=Float32)
    g=ForwardDiff.gradient(x -> parameter_energy(s,p,(;amplitudes=x);neighbourlist=w),theta)
    @test g ≈ oracle_gradient rtol=2e-5
    @test eltype(w.x)===Float32
end

@testset "Lennard-Jones parameter derivatives" begin
    p,_=ArgonComparison.argon_models()
    params=AC.get_parameters(p)
    theta=[params.epsilon,params.rmin]
    s=ArgonComparison.argon_pair(4.2)
    f=x -> ustrip(u"eV",AC.calculate(AC.Energy(),s,p,(epsilon=x[1],rmin=x[2])).energy)
    r=4.2; rc=p.radius; t=p.rmin/r; tc=p.rmin/rc
    df(x)=12*(x^11-x^5)
    ddf(x)=132*x^10-60*x^4
    exact=[f(theta)/p.epsilon,p.epsilon*((df(t)-df(tc))/r-(t-tc)*ddf(tc)/rc)]
    @test ForwardDiff.gradient(f,theta) ≈ exact rtol=1e-11
    @test all(isfinite,ForwardDiff.hessian(f,theta))
    @test AC.calculate(AC.Energy(),s,p;parameters=params).energy ≈ AC.potential_energy(s,p)
    g=energy_parameter_gradient(s,p,params)
    @test [g.epsilon,g.rmin] ≈ exact
    for sys in (ArgonComparison.argon_fcc(4),ArgonComparison.argon_pair(9.))
        f=x -> ustrip(AC.calculate(AC.Energy(),sys,p;
            parameters=(epsilon=x[1],rmin=x[2]),width=8,parallel=true).energy)
        ad=ForwardDiff.gradient(f,theta)
        fd=map(eachindex(theta)) do k
            h=1e-5*theta[k]; plus=copy(theta); minus=copy(theta)
            plus[k]+=h; minus[k]-=h
            (f(plus)-f(minus))/(2h)
        end
        @test ad ≈ fd atol=1e-8 rtol=2e-7
        @test Zygote.gradient(f,theta)[1] ≈ ad
    end
end

@testset "ChainRules energy pullback and calculate integration" begin
    s=fixture(19;cluster=true); p=potential(); params=AC.get_parameters(p)
    w=Workspace(width=8)
    expected=energy_parameter_gradient(s,p,params).amplitudes
    y,pullback=ChainRulesCore.rrule(parameter_energy,s,p,params;neighbourlist=w)
    @test y ≈ parameter_energy(s,p,params)
    tangent=pullback(2.)
    @test tangent[1] isa NoTangent
    @test tangent[2] isa ChainRulesCore.NotImplemented
    @test tangent[3] isa ChainRulesCore.NotImplemented
    @test tangent[4].amplitudes ≈ 2expected
    @test pullback(ZeroTangent())[4] isa ZeroTangent
    @test pullback(ChainRulesCore.Thunk(() -> 3.))[4].amplitudes ≈ 3expected
    # Pullbacks retain derivatives even after caller-owned mutable inputs change.
    saved=copy(params.amplitudes)
    params.amplitudes .*= 2
    s.xyz[1,1]+=0.3
    AC.calculate(AC.Energy(),s,p;neighbourlist=w)
    @test pullback(1.)[4].amplitudes ≈ expected
    params.amplitudes .= saved
    for positional in (false,true)
        f=x -> begin
            ps=(;amplitudes=x)
            r=positional ? AC.calculate(AC.Energy(),s,p,ps) :
                           AC.calculate(AC.Energy(),s,p;parameters=ps)
            ustrip(u"eV",r.energy)
        end
        @test Zygote.gradient(f,saved)[1] ≈ ForwardDiff.gradient(f,saved)
    end
end


@testset "multi-image parameter derivatives" begin
    H=[1.7 0.6 0.2;0 1.8 0.4;0 0 1.6]
    s=TestSystem(H*[0.1 0.4 0.8;0.2 0.6 0.3;0.1 0.7 0.5],H,(true,true,true),[:H,:He,:H])
    p=potential(); theta=AC.get_parameters(p).amplitudes
    expected=[image_oracle(s,with_parameters(p,(;amplitudes=Float64.(1:3 .== k)));shell=3).energy for k=1:3]
    for threshold in (0,32), width in (1,4), parallel in (false,true)
        f=x -> parameter_energy(s,p,(;amplitudes=x);width,parallel,bruteforce_threshold=threshold)
        @test ForwardDiff.gradient(f,theta) ≈ expected
        @test Zygote.gradient(f,theta)[1] ≈ expected
    end
end
