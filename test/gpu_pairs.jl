# Standalone: julia --project=<GPU environment> test/gpu_pairs.jl [--amdgpu]
# Without --amdgpu this explicitly tests KA CPU kernels, not a GPU.
using PotentialFrameWork, Test, Unitful, AtomsBase
import AtomsCalculators as AC
isdefined(@__MODULE__,:TestSystem) || include("fixtures.jl")
isdefined(@__MODULE__,:NL) || @eval const NL=PotentialFrameWork.NeighbourList
isdefined(@__MODULE__,:PP) || @eval const PP=PotentialFrameWork.PairPotentials
const gpu_test_backend=if "--amdgpu" in ARGS
    @eval using AMDGPU
    AMDGPU.functional() || error("AMDGPU was requested but is unusable")
    dev=only(filter(d->occursin("6800 XT",string(d)),AMDGPU.devices()))
    AMDGPU.device!(dev)
    AMDGPU.versioninfo()
    NL.KA.get_backend(AMDGPU.ROCArray(zeros(1)))
else
    @info "KA CPU kernel validation only"
    NL.KA.CPU()
end
function compare(s,p,w;parameters=nothing,request=(AC.Energy(),AC.Forces(),AC.Virial()))
    y=PP.evaluate_gpu!(w,s,p;parameters,request)
    T=eltype(w.x); tol=T===Float32 ? 3e-5 : 3e-11
    ref=evaluate(s,p;parameters,request,precision=T,width=1)
    for name in keys(y)
        @test ustrip.(getproperty(y,name)) ≈ ustrip.(getproperty(ref,name)) atol=tol rtol=tol
    end
    if p isa SmoothPair && !isempty(request)
        model=isnothing(parameters) ? p : with_parameters(p,parameters)
        converted=TestSystem(Float64.(T.(s.xyz)),Float64.(T.(s.H)),s.pbc,s.labels)
        converted_model=SmoothPair(model.labels,Float64.(T.(model.amplitudes)),Float64(T(model.radius)))
        independent=oracle(converted,converted_model)
        for name in keys(y)
            @test ustrip.(getproperty(y,name)) ≈ getproperty(independent,name) atol=tol rtol=tol
        end
    end
    @test eltype(w.x)===T && eltype(w.force)===T && eltype(w.partials)===T
    y
end
@testset "GPU radial pairs ($(typeof(gpu_test_backend)))" begin
 for T in (Float32,Float64), strategy in (:central,:candidate)
  w=NL.GPUWorkspace(gpu_test_backend;precision=T,config=NL.GPUConfig(;strategy))
  for n in (0,1,2,7,65,257)
   for pb in ((false,false,false),(true,false,true),(true,true,true))
    s=fixture(n;pbc=pb,cluster=true)
    compare(s,potential(),w)
   end
  end
  s=fixture(23)
  for mask=0:7
   props=Tuple(p for (bit,p) in enumerate((AC.Energy(),AC.Forces(),AC.Virial())) if (mask>>(bit-1))&1==1)
   compare(s,potential(),w;request=props)
  end
  y=compare(s,potential(),w); saved=copy(y.forces)
  ps=AC.get_parameters(potential()); ps.amplitudes .*= 2
  y2=compare(s,potential(),w;parameters=ps)
  @test y2.energy≈2*y.energy
  @test y.forces==saved
  perm=reverse(1:length(s)); sp=TestSystem(s.xyz[:,perm],s.H,s.pbc,s.labels[perm])
  yp=compare(sp,potential(),w)
  @test yp.forces≈y.forces[:,perm] rtol=(T===Float32 ? 3e-5 : 3e-11)
  # A two-cell periodic axis must visit each wrapped cell once.
  s=fixture(65;L=4.1); compare(s,potential(),w)
  s.xyz .+= [20.5,-12.3,8.2]; compare(s,potential(),w)
  # Named/isotopic identities use all 64 bits, including strict isequal identity.
  labels=[ChemicalSpecies(:C),ChemicalSpecies(:C13),ChemicalSpecies(:C;atom_name=:CA)]
  for label in labels; @test isequal(reinterpret(ChemicalSpecies,reinterpret(UInt64,label)),label); end
  p=SmoothPair(labels,[1. .4 .7;.4 2. .2;.7 .2 3.],2.)
  s=TestSystem([0. .5 1.;0. 0. 0.;0. 0. 0.],zeros(3,3),(false,false,false),labels)
  compare(s,p,w)
  s.xyz .= 0; compare(s,p,w) # regular distinct coincident origin
  @test_throws ArgumentError PP.evaluate_gpu!(w,s,potential())
  bad=fixture(2;L=4.); @test_throws ArgumentError PP.evaluate_gpu!(w,bad,potential())
  bad.H[1,2]=1; @test_throws ArgumentError PP.evaluate_gpu!(w,bad,potential())
  sparse=fixture(7;pbc=(false,false,false)); sparse.xyz[:,end].+=1e6
  compare(sparse,potential(),w)
  # Native input Float32 -> Float64 and Float64 -> Float32.
  s=fixture(7); compare(TestSystem(Float32.(s.xyz),Float32.(s.H),s.pbc,s.labels),potential(),w)
  lj=PP.lennard_jones(;sigma=.8,radius=2.)
  s=TestSystem([0. 1.2;0. .2;0. 0.],zeros(3,3),(false,false,false),[:Ar,:Ar])
  compare(s,lj,w)
  s.xyz .= 0
  @test_throws DomainError PP.evaluate_gpu!(w,s,lj)
  @test_throws DomainError PP.evaluate_gpu!(w,s,lj;request=(AC.Forces(),))
  @test_throws DomainError PP.evaluate_gpu!(w,s,lj;request=(AC.Virial(),))
  singular=PP.prepare_gpu_execution(w,s,lj)
  @test_throws DomainError PP.execute_gpu!(singular;request=(AC.Forces(),))
 end
end
@testset "Skin trajectories and invalidation" begin
 for T in (Float32,Float64),strategy in (:central,:candidate)
  w=NL.GPUWorkspace(gpu_test_backend;precision=T,config=NL.GPUConfig(;strategy,skin=.4))
  approaching=TestSystem([0. 2.2;0. 0.;0. 0.],zeros(3,3),(false,false,false),[:H,:He])
  compare(approaching,potential(),w); initial=w.rebuilds
  approaching.xyz[1,1]+=.15; approaching.xyz[1,2]-=.15
  pair_result=compare(approaching,potential(),w)
  @test pair_result.energy>0u"eV"
  @test w.rebuilds==initial
  approaching.xyz[1,1]+=.05
  compare(approaching,potential(),w)
  @test w.rebuilds>initial
  s=fixture(65)
  compare(s,potential(),w); base=w.rebuilds
  for dx in (.01,.08,.199,.2,.21,.3,-.01)
   s.xyz[1,1]+=dx
   compare(s,potential(),w)
  end
  @test w.rebuilds>base
  base=w.rebuilds
  s.xyz[1,1]+=s.H[1,1] # ambiguous wrapped jump conservatively rebuilds
  compare(s,potential(),w); @test w.rebuilds>base
  s.H[1,1]+=.1; compare(s,potential(),w)
  s.labels.=reverse(s.labels); compare(s,potential(),w)
  s=TestSystem(s.xyz[:,end:-1:1],s.H,s.pbc,s.labels[end:-1:1]); compare(s,potential(),w)
  compare(s,SmoothPair([:H,:He],[1. .7;.7 1.3],1.8),w)
 end
end
@testset "Finite difference force/strain" begin
 s=fixture(7;pbc=(false,false,false)); p=potential()
 w=NL.GPUWorkspace(gpu_test_backend)
 y=PP.evaluate_gpu!(w,s,p); h=1e-5
 s.xyz[1,2]+=h; ep=ustrip(PP.evaluate_gpu!(w,s,p).energy)
 s.xyz[1,2]-=2h; em=ustrip(PP.evaluate_gpu!(w,s,p).energy)
 s.xyz[1,2]+=h
 @test -(ep-em)/(2h)≈ustrip(y.forces[1,2]) atol=1e-8
 xp=copy(s.xyz); s.xyz[1,:].*=1+h; ep=ustrip(PP.evaluate_gpu!(w,s,p).energy)
 s.xyz.=xp; s.xyz[1,:].*=1-h; em=ustrip(PP.evaluate_gpu!(w,s,p).energy)
 @test -(ep-em)/(2h)≈ustrip(y.virial[1,1]) atol=1e-8
end

struct DeviceQuadratic{T}
    amplitude::T
end
function (p::DeviceQuadratic)(r,central,neighbors,request)
    request isa EnergyOnly && return (;energy=p.amplitude*r*r)
    (;energy=p.amplitude*r*r,derivative=2*p.amplitude*r)
end
PP.prepare_gpu_radial(p::DeviceQuadratic,::Type{T},ps) where T =
    DeviceQuadratic(PP.gpu_number(T,isnothing(ps) ? p.amplitude : ps.amplitude))
@testset "Custom radial device protocol and owning device results" begin
 s=fixture(13); p=RadialPotential(DeviceQuadratic(0.3),2.)
 for strategy in (:central,:candidate)
  w=NL.GPUWorkspace(gpu_test_backend;config=NL.GPUConfig(;strategy))
  y=compare(s,p,w)
  input=PP.prepare_gpu_execution(w,s,p)
  d=PP.execute_gpu!(input); NL.KA.synchronize(gpu_test_backend)
  saved=Array(d.forces)
  d2=PP.execute_gpu!(input;parameters=(;amplitude=.6)); NL.KA.synchronize(gpu_test_backend)
  @test Array(d2.forces)≈2*saved
  @test Array(d.forces)==saved
  @test only(Array(d.energy))*u"eV"≈y.energy
 end
end

@testset "GPU launch choices, request normalization and lifetime" begin
 s=fixture(129;cluster=true); p=potential()
 for workgroup in (32,128),tile in (1,2)
  w=NL.GPUWorkspace(gpu_test_backend;config=NL.GPUConfig(;strategy=:candidate,workgroup,tile))
  for _=1:3
   compare(s,p,w;request=(AC.Forces(),AC.Energy(),AC.Forces(),AC.Virial()))
  end
 end
 input=PP.prepare_gpu_execution(s,p,gpu_test_backend)
 saved=PP.execute_gpu!(input)
 PP.prepare_gpu_execution(input.workspace,s,p)
 @test_throws ArgumentError PP.execute_gpu!(input)
 @test_throws ArgumentError PP.evaluate_gpu!(input.workspace,s,p;request=(:unsupported,))
 @test_throws ArgumentError NL.GPUWorkspace(gpu_test_backend;precision=Float64,
     config=NL.GPUConfig(;strategy=:candidate,workgroup=256,tile=4))
 if Base.get_extension(PotentialFrameWork,:PotentialFrameWorkChainRulesCoreExt)!==nothing
  cr=Base.get_extension(PotentialFrameWork,:PotentialFrameWorkChainRulesCoreExt).ChainRulesCore
  @test_throws ArgumentError cr.rrule(PP.evaluate_gpu,s,p,gpu_test_backend)
  @test_throws ArgumentError cr.rrule(PP.execute_gpu!,input)
 end
end

# Exercise species accessors and strict ChemicalSpecies equality on device, not
# just the host-built parameter channel map.
struct SpeciesQuadratic{T}
    scale::T
end
function (p::SpeciesQuadratic)(r,center,neighbors,request)
    T=typeof(p.scale)
    a=p.scale*(isequal(center,neighbors[1]) ? one(T) : T(2))
    request isa EnergyOnly && return (;energy=a*r*r)
    (;energy=a*r*r,derivative=T(2)*a*r)
end
PP.prepare_gpu_radial(p::SpeciesQuadratic,::Type{T},ps) where T = SpeciesQuadratic(T(p.scale))
@testset "Device species accessor identity" begin
 labels=[ChemicalSpecies(:C),ChemicalSpecies(:C13),ChemicalSpecies(:C;atom_name=:CA)]
 s=TestSystem([0. .2 .4 .6 .8 1.;zeros(2,6)],zeros(3,3),(false,false,false),repeat(labels,2))
 for T in (Float32,Float64),strategy in (:central,:candidate)
  w=NL.GPUWorkspace(gpu_test_backend;precision=T,config=NL.GPUConfig(;strategy))
  compare(s,RadialPotential(SpeciesQuadratic(.1),2.),w)
 end
end
