import AtomsCalculators as AC
using StaticArrays: SVector, SMatrix

force_matrix(f)=[f[i][a] for a=1:3,i in eachindex(f)]
forces_agree(a,b)=isapprox(ustrip.(u"eV/Å",force_matrix(a)),ustrip.(u"eV/Å",force_matrix(b)))

struct EnergyOnlyPair <: AbstractPairPotential end
PotentialFrameWork.cutoff(::EnergyOnlyPair)=2.0
(::EnergyOnlyPair)(b::PairBlock,::EnergyOnly)=(;energy=zero(b.dx))

@testset "AtomsCalculators conformance" begin
    s=fixture(31;cluster=true)
    p=potential()
    AC.Testing.test_energy_forces_virial(s,p)
    # Radial adapters inherit the same interface, including all generated helpers.
    AC.Testing.test_energy_forces_virial(s,RadialPotential(SpeciesRadial(),2.))
end

@testset "AtomsCalculators numerical boundary and neighbourlist lifetime" begin
    p=potential()
    for n in (0,1,17,93), T in (Float64,Float32), parallel in (false,true)
        src=fixture(n;L=4.1,pbc=(true,false,true),cluster=true)
        s=TestSystem(T.(src.xyz),T.(src.H),src.pbc,src.labels)
        w=Workspace(width=4,precision=T)
        expected=evaluate(s,p;precision=T)
        r=AC.energy_forces_virial(s,p;neighbourlist=w,parallel)
        @test r.forces isa Vector{<:SVector{3}}
        @test r.virial isa SMatrix{3,3}
        @test typeof(ustrip(r.energy))===T
        @test r.energy ≈ expected.energy
        @test force_matrix(r.forces) ≈ expected.forces
        @test r.virial ≈ expected.virial
        # Utilities describe default output precision; an explicit Float32
        # workspace also works for empty systems (whose default is Float64).
        @test typeof(AC.forces(s,p))===typeof(AC.zero_forces(s,p))
        @test typeof(AC.potential_energy(s,p))===typeof(AC.zero_energy(s,p))
        @test typeof(AC.virial(s,p))===typeof(AC.zero_virial(s,p))
        @test eltype(AC.forces(s,p))===AC.promote_force_type(s,p)
        @test validate_schedule(w)
        @test AC.potential_energy(s,p;parallel) ≈ expected.energy
        @test force_matrix(AC.forces(s,p;parallel)) ≈ expected.forces
        @test AC.virial(s,p;parallel) ≈ expected.virial
        version=w.version
        ef=AC.energy_forces(s,p;neighbourlist=w,parallel)
        @test w.version==version+1 # one complete evaluation, not one per property
        @test ef.energy ≈ expected.energy
        @test force_matrix(ef.forces) ≈ expected.forces
        @test AC.energy_unit(p)==u"eV"
        @test AC.length_unit(p)==u"Å"
        @test AC.force_unit(p)==u"eV/Å"
    end
    s=fixture(37;L=4.1)
    w=Workspace(width=8,precision=Float32)
    for k=1:4
        s.xyz[1,1]+=0.1
        s.labels[2]=ChemicalSpecies(isodd(k) ? :H : :He)
        s.H[3,3]=4.1+0.1*k
        r=AC.energy_forces_virial(s,p;neighbourlist=w,parallel=true)
        agrees((;r.energy,forces=force_matrix(r.forces),r.virial),oracle(s,p);tol=3e-4)
    end
    @test_throws ArgumentError AC.forces(s,p;neighbourlist=:invalid)
    @test_throws ArgumentError AC.forces(s,p;neighbourlist=w,width=4)
    @test_throws ArgumentError AC.forces(s,p;neighbourlist=w,precision=Float64)
    @test AC.potential_energy(s,p;neighbourlist=w,width=8,precision=Float32) ≈
        AC.potential_energy(s,p;width=8,precision=Float32)
end

@testset "AtomsCalculators additive destinations and fused calculate" begin
    s=fixture(43;cluster=true); p=potential(); w=Workspace()
    expected=AC.energy_forces_virial(s,p)
    f=AC.zero_forces(s,p)
    @test AC.forces!(f,s,p;neighbourlist=w)===f
    @test forces_agree(f,expected.forces)
    @test AC.forces!(f,s,p;neighbourlist=w)===f
    @test forces_agree(f,2*expected.forces)
    version=w.version
    ef=AC.energy_forces!(f,s,p;neighbourlist=w)
    @test w.version==version+1
    @test ef.forces===f
    @test ef.energy ≈ expected.energy
    @test forces_agree(f,3*expected.forces)
    evf=AC.energy_forces_virial!(f,s,p;neighbourlist=w,parallel=true)
    @test evf.forces===f
    @test evf.virial ≈ expected.virial
    @test forces_agree(f,4*expected.forces)
    @test_throws DimensionMismatch AC.forces!(f[2:end],s,p)
    @test_throws DimensionMismatch AC.energy_forces!(f[2:end],s,p)
    @test_throws DimensionMismatch AC.energy_forces_virial!(f[2:end],s,p)
    props=(AC.Energy(),AC.Forces(),AC.Virial())
    state=(token=17,)
    version=w.version
    result=AC.calculate(props,s,p;state,neighbourlist=w,parallel=true)
    @test w.version==version+1
    @test keys(result)==(:energy,:forces,:virial,:state)
    @test result.energy ≈ expected.energy
    @test forces_agree(result.forces,expected.forces)
    @test result.virial ≈ expected.virial
    @test result.state===state
    positional=AC.calculate(props,s,p,nothing,state;neighbourlist=w)
    for other in (positional,AC.calculate(props,s,p,nothing;state,neighbourlist=w))
        @test other.energy ≈ result.energy
        @test forces_agree(other.forces,result.forces)
        @test other.virial ≈ result.virial
        @test other.state===state
    end
    for prop in props
        r=AC.calculate(prop,s,p;neighbourlist=w)
        key=first(keys(r))
        @test key===:forces ? forces_agree(r[key],expected[key]) : isapprox(r[key],expected[key])
        @test r.state===nothing
    end
    @test keys(AC.calculate((AC.Forces(),AC.Energy()),s,p))==(:forces,:energy,:state)
    @test AC.calculate((),s,p;state)==(;state)
    @test_throws ArgumentError AC.calculate(props,s,p;parameters=(radius=3.,))
    @test_throws ArgumentError AC.calculate((AC.Energy(),:unknown),s,p)
    @test AC.calculate(AC.Energy(),s,EnergyOnlyPair()).energy==0u"eV"
    @test AC.calculate((AC.Energy(),AC.Energy()),s,EnergyOnlyPair()).energy==0u"eV"
    @test AC.potential_energy(s,EnergyOnlyPair())==0u"eV"
end
