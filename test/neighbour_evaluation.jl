# Protocol conformance: callbacks below know only NI accessors. The independent
# O(N²) backend supports open cells only, intentionally uses different storage,
# and neither calls native neighbour search nor reads native view fields.
import AtomsCalculators as AC
using LinearAlgebra, StaticArrays
const NI=PotentialFrameWork.NeighbourInterface
struct ProtocolPair end
AC.length_unit(::ProtocolPair)=u"Å"
AC.energy_unit(::ProtocolPair)=u"eV"
NI.evaluation_kind(::ProtocolPair)=:pair
NI.search_cutoff(::ProtocolPair)=2.7
NI.species_channels(::ProtocolPair)=nothing
NI.prepare_potential(p::ProtocolPair,species;kwargs...)=p
NI.make_scratch(::ProtocolPair,n,T)=nothing
function NI.evaluate_pair!(out,p::ProtocolPair,b,ps,scratch,props)
    r=NI.distance(b)
    a=isnothing(ps) ? 0.8 : ps.a
    q=1-r*r/2.7^2
    any(p -> p isa AC.Energy,props) && NI.set_energy!(out,a*q*q*q)
    if any(p -> p isa Union{AC.Forces,AC.Virial},props)
        NI.set_radial_derivative!(out,-6a*q*q*r/2.7^2)
    end
    nothing
end
struct AnalyticSite end
AC.length_unit(::AnalyticSite)=u"Å"
AC.energy_unit(::AnalyticSite)=u"eV"
NI.evaluation_kind(::AnalyticSite)=:site
NI.search_cutoff(::AnalyticSite)=2.7
NI.species_channels(::AnalyticSite)=ChemicalSpecies.([:He,:H]) # deliberately reversed
NI.prepare_potential(p::AnalyticSite,species;kwargs...)=p
NI.make_scratch(::AnalyticSite,n,T)=nothing
# True angular many-body term: squared dot products of distinct edges.
function NI.evaluate_site!(out,::AnalyticSite,site,ps,scratch,props)
    energy=0.0
    for j=1:NI.neighbour_count(site)
        NI.active(site,j) || continue
        d=NI.displacement(site,j); g=zero(d)
        for k=1:NI.neighbour_count(site)
            (j==k || !NI.active(site,k)) && continue
            v=NI.displacement(site,k); c=dot(d,v)
            energy+=ps.a*c*c/2
            g+=2ps.a*c*v
        end
        any(p -> p isa Union{AC.Forces,AC.Virial},props) && NI.set_gradient!(out,j,g)
    end
    any(p -> p isa AC.Energy,props) && NI.set_energy!(out,energy)
    nothing
end
struct OracleSite
    atom::Int
    label::ChemicalSpecies
    edges::Vector{Tuple{ChemicalSpecies,SVector{3,Float64}}}
end
NI.central_index(s::OracleSite)=s.atom
NI.central_species(s::OracleSite)=s.label
NI.neighbour_count(s::OracleSite)=length(s.edges)
NI.active(s::OracleSite,k)=true
NI.displacement(s::OracleSite,k)=last(s.edges[k])
NI.neighbour_species(s::OracleSite,k)=first(s.edges[k])
struct OraclePair
    delta::NTuple{3,Float64}
end
NI.displacement(b::OraclePair)=b.delta
NI.distance(b::OraclePair)=sqrt(sum(abs2,b.delta))
mutable struct OracleOutput
    value::Float64
    gradients::Dict{Int,SVector{3,Float64}}
    pair_gradient::NTuple{3,Float64}
    radial_derivative::Float64
end
OracleOutput(v,g,p)=OracleOutput(v,g,p,NaN)
NI.set_radial_derivative!(o::OracleOutput,g)=(o.radial_derivative=g;nothing)
NI.set_energy!(o::OracleOutput,v)=(o.value=v;nothing)
NI.set_gradient!(o::OracleOutput,j,g)=(o.gradients[j]=g;nothing)
NI.set_gradient!(o::OracleOutput,g)=(o.pair_gradient=g;nothing)
function oracle_protocol(sys,p,ps;centers=1:length(sys))
    any(sys.pbc) && throw(ArgumentError("reference backend supports open geometry only"))
    E=0.; F=zeros(3,length(sys)); W=zeros(3,3)
    for i in centers
        js=[j for j=1:length(sys) if j!=i && norm(sys.xyz[:,j]-sys.xyz[:,i])<NI.search_cutoff(p)]
        out=OracleOutput(NaN,Dict{Int,SVector{3,Float64}}(),(NaN,NaN,NaN))
        if NI.evaluation_kind(p)===:site
            site=OracleSite(i,sys.labels[i],[(sys.labels[j],SVector{3}(sys.xyz[:,j]-sys.xyz[:,i])) for j in js])
            NI.evaluate_site!(out,p,site,ps,nothing,(AC.Energy(),AC.Forces(),AC.Virial()))
            E+=out.value
            for (k,j) in enumerate(js)
                g=out.gradients[k]; d=NI.displacement(site,k)
                F[:,i]+=g; F[:,j]-=g; W-=d*transpose(g)
            end
        else
            for j in js
                j>i || continue
                d=SVector{3}(sys.xyz[:,j]-sys.xyz[:,i])
                NI.evaluate_pair!(out,p,OraclePair(Tuple(d)),ps,nothing,(AC.Energy(),AC.Forces(),AC.Virial()))
                g=isnan(out.radial_derivative) ? SVector(out.pair_gradient) : out.radial_derivative*d/norm(d)
                E+=out.value
                F[:,i]+=g; F[:,j]-=g; W-=d*transpose(g)
            end
        end
    end
    (;energy=E,forces=F,virial=W)
end
struct BatchProbe
    direct::Base.RefValue{Int}
    batches::Base.RefValue{Int}
    order::Vector{Int}
end
AC.length_unit(::BatchProbe)=u"Å"
AC.energy_unit(::BatchProbe)=u"eV"
NI.evaluation_kind(::BatchProbe)=:site
NI.search_cutoff(::BatchProbe)=2.7
NI.species_channels(::BatchProbe)=NI.species_channels(AnalyticSite())
NI.prepare_potential(p::BatchProbe,species;kwargs...)=p
NI.make_scratch(::BatchProbe,n,T)=nothing
function NI.evaluate_site!(out,p::BatchProbe,site,ps,scratch,props)
    p.direct[]+=1; push!(p.order,NI.central_index(site))
    NI.evaluate_site!(out,AnalyticSite(),site,ps,scratch,props)
end
function NI.evaluate_sites!(out,p::BatchProbe,batch,ps,scratch,props)
    p.batches[]+=1
    for k=1:NI.site_count(batch)
        site=NI.site_at(batch,k); push!(p.order,NI.central_index(site))
        NI.evaluate_site!(NI.output_at(out,k),AnalyticSite(),site,ps,scratch,props)
    end
    nothing
end

@testset "experimental neighborhood protocol: independent backends and requests" begin
    sys=TestSystem([0. 0.9 1.4 4.9 5.2;0.1 0.4 1.1 0.2 0.7;0.2 0.8 0.3 0.5 1.1],
        12Matrix{Float64}(I,3,3),(false,false,false),[:H,:He,:H,:He,:H])
    ps=(;a=0.8)
    props=(AC.Energy(),AC.Forces(),AC.Virial())
    for bits=1:7
        request=Tuple(props[k] for k=1:3 if !iszero(bits & (1<<(k-1))))
        pairref=oracle_protocol(sys,ProtocolPair(),ps)
        result=evaluate(sys,ProtocolPair();request,parameters=ps)
        for (prop,key,unit) in zip(props,(:energy,:forces,:virial),(u"eV",u"eV/Å",u"eV"))
            if prop in request
                @test ustrip.(unit,result[key]) ≈ pairref[key]
            else
                @test !haskey(result,key)
            end
        end
        for width in (1,4), batch_size in (1,2,4,7), parallel in (false,true)
            input=NL.prepare_site_execution(AnalyticSite(),sys;width,batch_size,slots=3)
            w=input.execution
            for slot in w.slots, out in slot.outputs
                fill!(out.gradients,SVector(NaN,NaN,NaN)); out.energy=NaN
            end
            NL.execute!(w,input.positions,input.cell,ps;properties=request,parallel)
            result=NL.site_result(w); ref=oracle_protocol(sys,AnalyticSite(),ps)
            for (prop,key) in zip(props,(:energy,:forces,:virial))
                @test result[key] ≈ (prop in request ? ref[key] : zero(ref[key])) atol=1e-11
            end
        end
    end
    @test evaluate(sys,ProtocolPair();request=())==NamedTuple()
    @test keys(evaluate(sys,ProtocolPair();request=AC.Forces()))==(:forces,)
    @test keys(evaluate(sys,ProtocolPair();request=AC.Energy()))==(:energy,)
    @test_throws ArgumentError evaluate(sys,ProtocolPair();request=(:energy,))
    @test evaluate(sys,ProtocolPair();request=(AC.Virial(),AC.Energy(),AC.Virial())) ==
        evaluate(sys,ProtocolPair();request=(AC.Energy(),AC.Virial()))
    @test_throws ArgumentError NL.prepare_site_execution(AnalyticSite(),sys;batch_size=0)
    @test_throws ArgumentError NL.prepare_site_execution(AnalyticSite(),sys;batch_size=-1)
    # Independent numerical derivatives of the angular callback, including strain.
    ref=oracle_protocol(sys,AnalyticSite(),ps); h=1e-6
    for i=1:length(sys),a=1:3
        xp=deepcopy(sys); xm=deepcopy(sys); xp.xyz[a,i]+=h; xm.xyz[a,i]-=h
        fd=(oracle_protocol(xp,AnalyticSite(),ps).energy-oracle_protocol(xm,AnalyticSite(),ps).energy)/(2h)
        @test fd ≈ -ref.forces[a,i] atol=1e-7
    end
    strain=[0. 1 0;0 0 0;0 0 0]
    xp=TestSystem((I+h*strain)*sys.xyz,sys.H,sys.pbc,sys.labels)
    xm=TestSystem((I-h*strain)*sys.xyz,sys.H,sys.pbc,sys.labels)
    @test (oracle_protocol(xp,AnalyticSite(),ps).energy-oracle_protocol(xm,AnalyticSite(),ps).energy)/(2h) ≈ -sum(ref.virial.*transpose(strain)) atol=1e-7
    for batch_size in (1,2,4)
        probe=BatchProbe(Ref(0),Ref(0),Int[])
        input=NL.prepare_site_execution(probe,sys;batch_size)
        w=input.execution
        for centers in (nothing,Int[],[3],[5,2,1])
            empty!(probe.order); probe.direct[]=0; probe.batches[]=0
            NL.execute!(w,input.positions,input.cell,ps;centers)
            expected=isnothing(centers) ? collect(1:length(sys)) : centers
            @test probe.order==expected
            @test probe.direct[]==(batch_size==1 ? length(expected) : 0)
            @test probe.batches[]==(batch_size==1 ? 0 : cld(length(expected),batch_size))
            @test NL.site_result(w).forces ≈ oracle_protocol(sys,AnalyticSite(),ps;centers=expected).forces
        end
        @test_throws ArgumentError NL.execute!(w,input.positions,input.cell,ps;centers=[1,1])
        @test_throws ArgumentError NL.execute!(w,input.positions,input.cell,ps;centers=[0])
        saved=NL.site_result(w)
        NL.execute!(w,input.positions,input.cell,(;a=1.6))
        @test NL.site_result(w).energy ≈ 2oracle_protocol(sys,AnalyticSite(),ps).energy
        @test saved.forces !== w.force
        w.environment.version+=UInt(1)
        @test_throws ArgumentError NL.execute!(w,input.positions,input.cell,ps)
    end
    periodic=TestSystem(sys.xyz,sys.H,(true,false,false),sys.labels)
    @test_throws ArgumentError oracle_protocol(periodic,AnalyticSite(),ps)
end

struct UnsupportedUnits end
AC.length_unit(::UnsupportedUnits)=u"nm"
AC.energy_unit(::UnsupportedUnits)=u"eV"
NI.evaluation_kind(::UnsupportedUnits)=:site

@testset "ACE protocol batches, periodic images, precision and compatibility" begin
    A=PotentialFrameWork.ACE
    model,ps,st=A.load_resolved(joinpath(@__DIR__,"../data/ace_ref/resolved_model_v1.toml"))
    H=[3.7 0.7 0.1;0 3.8 0.3;0 0 3.6]
    sys=TestSystem([0.1 2.3 1.1;0.2 0.7 2.4;0.3 1.2 1.5],H,(true,true,true),[:Ti,:Al,:Ti])
    @test_throws ArgumentError NL.prepare_site_execution(UnsupportedUnits(),sys)
    @test_throws ArgumentError NL.prepare_site_execution(model,sys;simd=false,fused=true)
    for T in (Float32,Float64),width in (1,4),batch_size in (1,2,5)
        input=A.ace_inputs(sys,model;environment=A.FullEnvironment(;width))
        w=A.CPUExecutor(model,input.environment;precision=T,batch_size,slots=3)
        for centers in (nothing,Int[],[2],[3,1]),parallel in (false,true)
            ref=A.ace_evaluate(input.positions,model,ps,input.environment,input.cell;centers)
            for bits=1:7
                allprops=(AC.Energy(),AC.Forces(),AC.Virial())
                props=Tuple(allprops[k] for k=1:3 if !iszero(bits & (1<<(k-1))))
                A.execute_prepared!(w,input.positions,input.cell,ps;centers,parallel,properties=props)
                y=A.prepared_result(w); tol=T===Float32 ? 4e-4 : 3e-10
                for (p,key) in zip(allprops,(:energy,:forces,:virial))
                    @test y[key] ≈ (p in props ? ref[key] : zero(ref[key])) atol=tol rtol=tol
                end
            end
        end
        # Loop fallback must preserve the angular backend and fusion options.
        for angular_backend in (:sphericart,:generated,:block),fused in (false,true)
            wb=A.CPUExecutor(model,input.environment;precision=T,batch_size,angular_backend,slots=3)
            A.execute_prepared!(wb,input.positions,input.cell,ps;fused)
            ref=A.ace_evaluate(input.positions,model,ps,input.environment,input.cell)
            @test A.prepared_result(wb).forces ≈ ref.forces atol=4e-4 rtol=4e-4
        end
        @test sizeof(w.force)==3length(sys)*sizeof(T)
    end
    for batch_size in (1,2,5)
        input=NL.prepare_site_execution(model,sys;batch_size,width=4)
        w=input.execution
        for scale in (1.,0.7)
            updated=A.install_readout(ps,scale*A.linear_parameters(ps))
            NL.execute!(w,input.positions,input.cell,updated)
            ref=A.ace_evaluate(input.positions,model,updated,w.environment,input.cell)
            @test NL.site_result(w).energy ≈ ref.energy
            @test NL.site_result(w).forces ≈ ref.forces atol=2e-10
        end
        calc=A.ACECalculator(model,ps;batch_size)
        for properties in ((AC.Forces(),),(AC.Virial(),),(AC.Virial(),AC.Energy(),AC.Energy()))
            result=AC.calculate(properties,sys,calc;state=(;marker=3))
            @test result.state==(;marker=3)
            @test keys(result)==(Tuple(p isa AC.Energy ? :energy : p isa AC.Forces ? :forces : :virial for p in unique(properties))...,:state)
        end
    end
    self=TestSystem(zeros(3,1),2.8Matrix{Float64}(I,3,3),(true,true,true),[:Ti])
    input=NL.prepare_site_execution(model,self;batch_size=2,width=4)
    NL.execute!(input.execution,input.positions,input.cell,ps)
    y=NL.site_result(input.execution)
    @test all(iszero,y.forces)
    @test !iszero(y.energy)
    @test !all(iszero,y.virial)
    @test_throws ArgumentError A.CPUExecutor(model,input.execution.environment;batch_size=0)
end

@testset "empty neighborhoods, layout lifetime and transformed batches" begin
    for n in (0,1), T in (Float32,Float64), width in (1,4), batch_size in (1,2,5)
        sys=TestSystem(zeros(3,n),10Matrix{Float64}(I,3,3),(false,false,false),fill(:H,n))
        input=NL.prepare_site_execution(AnalyticSite(),sys;precision=T,width,batch_size)
        NL.execute!(input.execution,input.positions,input.cell,(;a=0.8))
        y=NL.site_result(input.execution)
        @test iszero(y.energy)
        @test all(iszero,y.forces)
        @test all(iszero,y.virial)
        input.execution.environment.width=width==1 ? 4 : 1
        @test_throws ArgumentError NL.execute!(input.execution,input.positions,input.cell,(;a=0.8))
    end
    A=PotentialFrameWork.ACE
    m,ps,st=A.load_resolved(joinpath(@__DIR__,"../data/ace_ref/resolved_model_v1.toml"))
    selected=A.BasisSelection(m;indices=1:20)
    pruned=A.prune_basis(m,ps,selected)
    cases=((A.mask_basis(m,selected),ps),(pruned.model,pruned.parameters),(A.project_species(m,[:Ti],1:2),ps))
    sys=TestSystem([0.1 2.3 1.1;0.2 0.7 2.4;0.3 1.2 1.5],8Matrix{Float64}(I,3,3),(false,false,false),[:Ti,:Al,:Ti])
    for (model,parameters) in cases,batch_size in (2,5)
        input=A.ace_inputs(sys,model;environment=A.FullEnvironment(width=4))
        w=A.CPUExecutor(model,input.environment;batch_size)
        A.execute_prepared!(w,input.positions,input.cell,parameters;centers=[3,1])
        ref=A.ace_evaluate(input.positions,model,parameters,input.environment,input.cell;centers=[3,1])
        y=A.prepared_result(w)
        @test y.energy ≈ ref.energy
        @test y.forces ≈ ref.forces atol=2e-10
        @test y.virial ≈ ref.virial atol=2e-10
        input.environment.width=1
        @test_throws ArgumentError A.execute_prepared!(w,input.positions,input.cell,parameters)
    end
end

# A radial callback can be consumed by a backend that exposes no displacement.
struct DistanceOnlyPair
    radius::Float64
end
NI.distance(b::DistanceOnlyPair)=b.radius
@testset "distance-first radial pair outputs and optional Cartesian access" begin
    out=OracleOutput(NaN,Dict{Int,SVector{3,Float64}}(),(NaN,NaN,NaN))
    NI.evaluate_pair!(out,ProtocolPair(),DistanceOnlyPair(0.7),(;a=0.8),nothing,(AC.Energy(),AC.Forces()))
    q=1-0.7^2/2.7^2
    @test out.value ≈ 0.8q^3
    @test out.radial_derivative ≈ -6*0.8q^2*0.7/2.7^2
    @test all(isnan,out.pair_gradient) # radial callback supplied no orientation
    for T in (Float32,Float64),width in (1,4,8),parallel in (false,true)
        sys=fixture(17;cluster=true)
        # Nontrivial radial polynomial has a regular origin and analytic derivative.
        radial=(r,si,sj,req) -> req isa EnergyOnly ? (;energy=(1-r*r/4)^3) :
            (;energy=(1-r*r/4)^3,derivative=-1.5r*(1-r*r/4)^2)
        p=RadialPotential(radial,2.)
        legacy=PP.LegacyPair(p)
        w=Workspace(;width,precision=T)
        NL.prepare!(w,sys,cutoff(p))
        props=(AC.Energy(),AC.Forces(),AC.Virial())
        prepared=NI.prepare_potential(p,w.species_map)
        NL.execute!(w,legacy,props;parallel)
        reference=PP.outputs(w,EnergyForcesVirial())
        NL.execute!(w,prepared,props;parallel)
        y=PP.outputs(w,EnergyForcesVirial())
        tol=T===Float32 ? 3e-5 : 1e-12
        @test y.energy ≈ reference.energy rtol=tol
        @test y.forces ≈ reference.forces rtol=tol atol=tol*u"eV/Å"
        @test y.virial ≈ reference.virial rtol=tol atol=tol*u"eV"
        # New default native callbacks use radial outputs; old Cartesian callbacks
        # still work and are exercised by the anisotropic regression suite.
        smooth=NI.prepare_potential(potential(),w.species_map)
        @test smooth isa PP.PreparedSmooth
    end
    # Explicit parameters retain the model's original triangular coordinates,
    # even when provider species channels occur in a different order.
    sys=fixture(9;cluster=true)
    model=SmoothPair([:He,:H],[0.3 0.7;0.7 0.8],2.)
    w=Workspace(width=4); NL.prepare!(w,sys,cutoff(model))
    prepared=NI.prepare_potential(model,w.species_map)
    props=(AC.Energy(),AC.Forces(),AC.Virial())
    for scale in (0.6,1.4)
        ps=(;amplitudes=scale*AC.get_parameters(model).amplitudes)
        NL.execute!(w,prepared,props;parameters=ps)
        y=PP.outputs(w,EnergyForcesVirial())
        ref=evaluate(sys,with_parameters(model,ps))
        @test y.energy ≈ ref.energy
        @test y.forces ≈ ref.forces
        @test y.virial ≈ ref.virial
    end
    @test_throws ArgumentError NL.execute!(w,prepared,props;parameters=(;amplitudes=[1.]))
    sys=TestSystem(zeros(3,2),10Matrix{Float64}(I,3,3),(false,false,false),[:H,:H])
    radial=(r,si,sj,req) -> req isa EnergyOnly ? (;energy=r*r) : (;energy=r*r,derivative=2r)
    y=evaluate(sys,RadialPotential(radial,2.))
    @test iszero(y.energy)
    @test all(iszero,y.forces)
    @test all(iszero,y.virial)
end
