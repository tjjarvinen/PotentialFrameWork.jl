include(joinpath(@__DIR__,"..","benchmark","lennard_jones","common.jl"))
using .ArgonComparison

@testset "Lennard-Jones agreement with EmpiricalPotentials $(pkgversion(ArgonComparison.EP))" begin
    ours,reference=argon_models()
    # Analytic inverse-distance tilt vs upstream's ForwardDiff force, including
    # both sides of and the exact cutoff. No per-atom output reordering in tests.
    for r in (3.3,3.8,5.,ours.radius-1e-5,ours.radius,ours.radius+1e-5,12.)
        s=argon_pair(r)
        ref=AC.energy_forces_virial(s,reference)
        for width in (1,4,8), parallel in (false,true)
            actual=AC.energy_forces_virial(s,ours;width,parallel)
            @test isapprox(numeric(actual).energy,numeric(ref).energy;rtol=2e-11,atol=1e-12)
            @test isapprox(numeric(actual).forces,numeric(ref).forces;rtol=2e-11,atol=1e-12)
            @test isapprox(numeric(actual).virial,numeric(ref).virial;rtol=2e-11,atol=1e-12)
        end
    end
    for (m,pbc) in ((1,(false,false,false)),(2,(false,false,false)),
                    (4,(true,true,true)),(4,(true,false,true)),(5,(true,true,true)))
        s=argon_fcc(m;periodic=pbc)
        original=copy(AtomsBase.position(s,:))
        ref=AC.energy_forces_virial(s,reference)
        for width in (1,4,8), parallel in (false,true)
            w=Workspace(;width)
            r=AC.energy_forces_virial(s,ours;neighbourlist=w,parallel)
            @test assert_agreement(r,ref) isa NamedTuple
            @test AC.potential_energy(s,ours;neighbourlist=w,parallel) ≈ ref.energy rtol=2e-10
            @test validate_schedule(w)
            @test numeric(r).energy<0
            @test maximum(abs,sum(numeric(r).forces;dims=2)) < 1e-11
        end
        @test AtomsBase.position(s,:)==original
    end
    h=1e-5; r=3.8
    ep=AC.potential_energy(argon_pair(r+h),ours)
    em=AC.potential_energy(argon_pair(r-h),ours)
    f=AC.forces(argon_pair(r),ours)[1][1]
    @test ustrip(u"eV/Å",f) ≈ ustrip(u"eV",ep-em)/(2h) rtol=1e-8
    @test_throws DomainError AC.forces(argon_pair(0.),ours)
    bad=ArgonComparison.AtomsSystems.generic_system([:Ne],[SVector(0.,0.,0.)*u"Å"])
    @test_throws ArgumentError AC.potential_energy(bad,ours)
end
