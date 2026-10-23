using PotentialFrameWork, Test, Unitful, AtomsBase
include("fixtures.jl")
const NL=PotentialFrameWork.NeighbourList
const PP=PotentialFrameWork.PairPotentials

struct AnisotropicPair
    radius::Float64
end
PotentialFrameWork.cutoff(p::AnisotropicPair)=p.radius
function (p::AnisotropicPair)(b::PairBlock,req)
    x,y,z=b.dx,b.dy,b.dz
    q=1-(x*x+y*y+z*z)/p.radius^2
    h=x*x+2*y*y+0.3*x*y
    energy=h*q*q*q
    req isa EnergyOnly && return (;energy)
    g=-6*h*q*q/p.radius^2
    (;energy,fx=(2*x+0.3*y)*q*q*q+g*x,
        fy=(4*y+0.3*x)*q*q*q+g*y,fz=g*z)
end

function agrees(result,ref;tol=1e-11)
    @test ustrip(u"eV",result.energy) ≈ ref.energy atol=tol rtol=tol
    if haskey(result,:forces)
        @test ustrip.(u"eV/Å",result.forces) ≈ ref.forces atol=tol rtol=tol
        @test ustrip.(u"eV",result.virial) ≈ ref.virial atol=tol rtol=tol
    end
end

struct SpeciesRadial end
function (::SpeciesRadial)(r,central::ChemicalSpecies,neighbours::NTuple{W,ChemicalSpecies},req) where W
    req isa EnergyOnly && return (;energy=r*r)
    (;energy=r*r,derivative=2*r)
end

@testset "ChemicalSpecies storage, callbacks and strict identity" begin
    labels=[ChemicalSpecies(:C),ChemicalSpecies(:C13),ChemicalSpecies(:C;atom_name=:CA)]
    @test labels[1]==labels[2] # AtomsBase wildcard matching must not merge parameters.
    @test !isequal(labels[1],labels[2])
    s=TestSystem([0. 0.6 1.2;0. 0.2 0.1;0. 0. 0.],
        [8. 0 0;0 8 0;0 0 8],(true,true,true),labels)
    p=SmoothPair(labels,[1. 0.7 0.4;0.7 2. 0.3;0.4 0.3 3.],2.)
    @test eltype(p.labels)===ChemicalSpecies
    for width in (1,4,8)
        w=Workspace(;width)
        agrees(evaluate!(w,s,p;parallel=true),oracle(s,p))
        @test eltype(w.species)===ChemicalSpecies
        @test eltype(w.species_map)===ChemicalSpecies
        @test length(w.species_map)==3
        @test isequal(w.species,s.labels[w.p])
        _,block=NL.geometry(w,1,2,3)
        @test block.central_species isa ChemicalSpecies
        @test block.neighbour_species isa NTuple{width,ChemicalSpecies}
        perm=[3,1,2]
        sp=TestSystem(s.xyz[:,perm],s.H,s.pbc,s.labels[perm])
        agrees(evaluate!(w,sp,p),oracle(sp,p))
        radial=RadialPotential(SpeciesRadial(),2.)
        scalar=evaluate(s,radial;width=1)
        vector=evaluate(s,radial;width)
        @test vector.energy ≈ scalar.energy
        @test vector.forces ≈ scalar.forces
    end
    @test_throws ArgumentError evaluate(s,SmoothPair([ChemicalSpecies(:C)],ones(1,1),2.))
    @test_throws ArgumentError SmoothPair([:C,ChemicalSpecies(:C)],ones(2,2),2.)
end

@testset "serial scalar oracle, species, units and permutation" begin
    p=potential()
    for n in (0,1,2,7,31,93), pbc in ((false,false,false),(true,true,true),(true,false,true)), cluster in (false,true)
        s=fixture(n;pbc,cluster)
        ref=oracle(s,p)
        w=Workspace(width=1)
        result=evaluate!(w,s,p)
        agrees(result,ref)
        agrees(evaluate!(w,s,p;request=EnergyOnly()),ref)
        @test w.p[w.inverse]==collect(1:n)
        @test sum(ustrip.(u"eV/Å",result.forces);dims=2) ≈ zeros(3,1) atol=1e-12
        @test s.xyz==fixture(n;pbc,cluster).xyz
    end
    s=fixture(35; pbc=(false,false,false))
    r=evaluate(s,p;width=1)
    perm=reverse(1:length(s))
    sp=TestSystem(s.xyz[:,perm],s.H,s.pbc,s.labels[perm])
    rp=evaluate(sp,p;width=1)
    @test rp.energy ≈ r.energy
    @test rp.forces ≈ r.forces[:,perm]
    translated=TestSystem(s.xyz .+ [3.2,-9.1,4.7],s.H,s.pbc,s.labels)
    @test evaluate(translated,p;width=1).energy ≈ r.energy
    # Explicit nm -> Å conversion, actual AtomsBase FlexibleSystem interface.
    native=AtomsBase.atomic_system([:H=>[0.,0.,0.]*u"nm",:He=>[0.1,0.,0.]*u"nm"],
        ([1.,0.,0.]*u"nm",[0.,1.,0.]*u"nm",[0.,0.,1.]*u"nm"),(true,false,true))
    labels=AtomsBase.species(native,:)
    np=SmoothPair(labels,[1.0 0.7;0.7 1.3],2.)
    @test evaluate(native,np;width=1).energy ≈ (0.7*(1-1/4)^3)*u"eV"
    iso=AtomsBase.isolated_system([:H=>[0.,0.,0.]*u"Å",:H=>[1.,0.,0.]*u"Å"])
    @test evaluate(iso,SmoothPair([AtomsBase.species(iso,1)],ones(1,1),2.)).energy ≈ (1-1/4)^3*u"eV"
end

@testset "force and virial derivatives" begin
    s=fixture(14;pbc=(false,false,false),cluster=true); p=potential(); h=1e-5
    r=evaluate(s,p;width=1)
    for i=1:length(s), a=1:3
        xp=copy(s.xyz); xm=copy(s.xyz); xp[a,i]+=h; xm[a,i]-=h
        ep=evaluate(TestSystem(xp,s.H,s.pbc,s.labels),p;request=EnergyOnly(),width=1).energy
        em=evaluate(TestSystem(xm,s.H,s.pbc,s.labels),p;request=EnergyOnly(),width=1).energy
        @test ustrip(u"eV/Å",r.forces[a,i]) ≈ -ustrip(u"eV",ep-em)/(2h) atol=1e-8
    end
    for pbc in ((false,false,false),(true,true,true))
        s=fixture(20;pbc)
        r=evaluate(s,p;width=1)
        for strain in ([1. 0 0;0 1 0;0 0 1], [0. 1 0;0 0 0;0 0 0])
            I3=[1. 0 0;0 1 0;0 0 1]
            ap=I3+h*strain; am=I3-h*strain
            hp=pbc[1] ? ap*s.H : s.H; hm=pbc[1] ? am*s.H : s.H
            ep=evaluate(TestSystem(ap*s.xyz,hp,pbc,s.labels),p;request=EnergyOnly(),width=1).energy
            em=evaluate(TestSystem(am*s.xyz,hm,pbc,s.labels),p;request=EnergyOnly(),width=1).energy
            expected=-sum(ustrip.(u"eV",r.virial).*transpose(strain))
            @test ustrip(u"eV",ep-em)/(2h) ≈ expected atol=1e-7
        end
    end
end

@testset "full Cartesian callback and nonsymmetric virial convention" begin
    s=fixture(8;cluster=true,pbc=(false,false,false)); p=AnisotropicPair(2.)
    r=evaluate(s,p); h=1e-5
    for i=1:length(s), a=1:3
        xp=copy(s.xyz); xm=copy(s.xyz); xp[a,i]+=h; xm[a,i]-=h
        ep=evaluate(TestSystem(xp,s.H,s.pbc,s.labels),p;request=EnergyOnly()).energy
        em=evaluate(TestSystem(xm,s.H,s.pbc,s.labels),p;request=EnergyOnly()).energy
        @test ustrip(u"eV/Å",r.forces[a,i]) ≈ -ustrip(u"eV",ep-em)/(2h) atol=1e-8
    end
    @test !isapprox(r.virial,transpose(r.virial);atol=1e-8*u"eV")
    I3=[1. 0 0;0 1 0;0 0 1]
    for strain in ([0. 1 0;0 0 0;0 0 0],[0. 0 0;0 0 1;0 0 0],I3)
        ep=evaluate(TestSystem((I3+h*strain)*s.xyz,s.H,s.pbc,s.labels),p;request=EnergyOnly()).energy
        em=evaluate(TestSystem((I3-h*strain)*s.xyz,s.H,s.pbc,s.labels),p;request=EnergyOnly()).energy
        @test ustrip(u"eV",ep-em)/(2h) ≈ -sum(ustrip.(u"eV",r.virial).*transpose(strain)) atol=1e-7
    end
end

@testset "SIMD widths, precision, tails, exact cutoff and safe inactive lanes" begin
    p=potential()
    for T in (Float64,Float32), width in (1,2,4,8,16), n in (0,1,3,9,47)
        s=fixture(n;L=4.1,cluster=true)
        agrees(evaluate(s,p;width,precision=T),oracle(s,p);tol=T==Float32 ? 3e-4 : 1e-11)
    end
    s=TestSystem([0. 2. 4. 4.1;0 0 0 0;0 0 0 0], [10. 0 0;0 10 0;0 0 10],
        (false,false,false),fill(:H,4))
    agrees(evaluate(s,p;width=8),oracle(s,p))
    edge=TestSystem([0. 2.;0 0;0 0],s.H,s.pbc,[:H,:H])
    @test iszero(evaluate(edge,p).energy)
    @test all(iszero,evaluate(edge,p).forces)
    edge.xyz[:,2].=edge.xyz[:,1]
    @test evaluate(edge,p).energy == 1u"eV"
    @test all(iszero,evaluate(edge,p).forces)
    # All padded/inactive geometry is nonzero before division. Energy-only adapter
    # never requests derivatives; error if that branch were evaluated.
    radial=(r,si,sj,req)-> req isa EnergyOnly ? (;energy=1/r) : (;energy=1/r,derivative=-1/(r*r))
    adapter=RadialPotential(radial,2.)
    for n in (1,2,3,7)
        s=fixture(n;L=4.1)
        a=evaluate(s,adapter;width=1)
        for width in (2,4,8)
            b=evaluate(s,adapter;width)
            @test isfinite(ustrip(b.energy))
            @test b.energy ≈ a.energy
            @test b.forces ≈ a.forces
        end
        @test evaluate(s,adapter;request=EnergyOnly()).energy ≈ a.energy
    end
end

@testset "periodic ownership and actual write footprints" begin
    for L in (4.01,6.1,8.1,11.), pbc in ((true,true,true),(true,false,true),(false,false,false))
        s=fixture(77;L,pbc,cluster=true); p=potential(); w=Workspace(width=4)
        evaluate!(w,s,p)
        @test validate_schedule(w)
        counts=Dict{Tuple{Int,Int},Int}()
        for task in eachindex(w.ranges)
            actual=Set{Int}()
            for i in w.ranges[task], cell in w.adjacency[task]
                r=w.ranges[cell]
                for j=max(first(r),i+1):4:last(r)
                    js,b=NL.geometry(w,i,j,last(r))
                    for k=1:4
                        NL.lane(b.active,k) || continue
                        a,bb=minmax(w.p[i],w.p[js[k]])
                        counts[(a,bb)]=get(counts,(a,bb),0)+1
                        push!(actual,task); push!(actual,cell)
                    end
                end
            end
            @test issubset(actual,Set(w.footprints[task]))
            # Check cells from actual force destinations, not just loop metadata.
            f0=copy(w.force); fill!(w.force,0)
            PP.cell_task!(w,bind_species(p,w.species_map),EnergyForcesVirial(),task)
            destcells=Set(c for c in eachindex(w.ranges) if any(!iszero,w.force[:,w.ranges[c]]))
            @test issubset(destcells,Set(w.footprints[task]))
            w.force.=f0
        end
        @test Set(keys(counts))==oracle(s,p).pairs
        @test all(==(1),values(counts))
        w.schedule_version-=UInt(1)
        @test !validate_schedule(w)
        @test_throws ErrorException PP.execute!(w,bind_species(p,w.species_map),EnergyForcesVirial())
    end
    s=fixture(25;L=4.1); p=potential()
    # Far outside the primary cell, exact faces, and crossing periodic boundaries.
    s.xyz[:,1].=[0.,4.1,8.2]; s.xyz[:,2].=[4.09,0.01,0.03]
    s.xyz[:,3].+= [8.2,-12.3,4.1]
    agrees(evaluate(s,p),oracle(s,p))
    for width in (1,4,8)
        w=Workspace(;width)
        for repeat=1:6
            s.xyz[1,2]+=0.02
            s.labels[2]=ChemicalSpecies(isodd(repeat) ? :H : :He)
            s.H[3,3]=4.1+0.01*repeat
            agrees(evaluate!(w,s,p),oracle(s,p))
        end
    end
end

@testset "parallel one-array execution and reused outputs" begin
    p=potential()
    for n in (0,1,5,121,513), cluster in (false,true), width in (1,4)
        s=fixture(n;L=n>200 ? 20. : 6.1,cluster)
        ref=oracle(s,p); w=Workspace(;width)
        for repeat=1:5
            agrees(evaluate!(w,s,p;parallel=true),ref)
            agrees(evaluate!(w,s,p;parallel=true,request=EnergyOnly()),ref)
            @test sizeof(w.force)==3*n*sizeof(Float64)
            @test validate_schedule(w)
        end
    end
end

@testset "clear unsupported input errors" begin
    p=potential(); s=fixture(5); s.H[:,2]=s.H[:,1]
    @test_throws ArgumentError evaluate(s,p)
    s=fixture(5); s.H[1,1]=Inf
    @test_throws ArgumentError evaluate(s,p)
    s=fixture(5); s.xyz[1,1]=NaN
    @test_throws ArgumentError evaluate(s,p)
    @test_throws ArgumentError evaluate(fixture(3),SmoothPair([:H],ones(1,1),2.))
    @test_throws ArgumentError Workspace(width=3)
    @test_throws ArgumentError SmoothPair([:H,:He],[1. 2.;3. 1.],2.)
    @test_throws ArgumentError NL.length_value(Float64,1.0)
    @test_throws ArgumentError NL.length_value(Float64,1u"eV")
end
include("periodic_images.jl")
include("neighbour_interface.jl")
include("neighbour_benchmarks.jl")
include("calculators.jl")
include("lennard_jones.jl")
include("parameter_ad.jl")
include("numerical_pair_ad.jl")
include("ace.jl")
include("ace_cpu_ad.jl")
include("ace_angular.jl")
println("Validation process: Julia ",VERSION,", threads=",Threads.nthreads())
