# Portable arithmetic check; --amdgpu requires actual RX 6800 XT execution.
using PotentialFrameWork, Test, Unitful, AtomsBase, LinearAlgebra, StaticArrays
import AtomsCalculators as AC
isdefined(@__MODULE__,:TestSystem) || include("fixtures.jl")
const GA=PotentialFrameWork.ACE
const GN=PotentialFrameWork.NeighbourList
const array_constructor=if "--amdgpu" in ARGS
    @eval using AMDGPU
    AMDGPU.functional() || error("AMDGPU requested but unavailable")
    AMDGPU.device!(only(filter(d->occursin("6800 XT",string(d)),AMDGPU.devices())))
    AMDGPU.versioninfo()
    AMDGPU.ROCArray
else
    @info "KA CPU arithmetic validation, not AMD hardware validation"
    Array
end
function gpu_bundle(m,p,T)
    b=GA.adapt_model(m,p;precision=T,array_type=array_constructor)
    array_constructor===Array ? (;model=GA.GPUACEModel(b.model,GA.pack_model(b.model,Array),T),parameters=b.parameters) : b
end
function cpu_reference(sys,m,p;centers=nothing)
    T=eltype(m.reference); input=GA.ace_inputs(GA.precision_geometry(sys,T),m;precision=T)
    ex=GA.CPUExecutor(m,input.environment;precision=T,angular_backend=:block)
    GA.execute_prepared!(ex,T.(input.positions),T.(input.cell),p;centers)
    GA.prepared_result(ex)
end
struct ScalarParameterRadial end
function (::ScalarParameterRadial)(v,d,r,p,c)
    q=max(one(r)-r/c.cutoff,zero(r))
    for n in eachindex(v)
        v[n]=p.rate*q*q/(n+1)
        isnothing(d) || (d[n]=-2*p.rate*q/(c.cutoff*(n+1)))
    end
    nothing
end
@testset "GPU scalar ACE" begin
    m,p,_=GA.load_resolved(joinpath(@__DIR__,"../data/ace_ref/resolved_model_v1.toml"))
    sys=TestSystem([0. 2.6 1.1 8.;0. .2 2.8 7.;0. .3 .1 6.],Matrix{Float64}(I,3,3)*16,(true,true,true),[:Ti,:Al,:Ti,:Al])
    for T in (Float64,Float32)
        b=gpu_bundle(m,p,T); cpu=GA.adapt_model(m,p;precision=T)
        tol=T===Float64 ? 2e-9 : 4e-3
        ref=cpu_reference(sys,cpu.model,cpu.parameters)
        for strategy in (:staged,:baseline,:compact), batch in (1,3,10)
            ex=GA.GPUExecutor(sys,b.model,b.parameters;strategy,batch_size=batch)
            for subset=0:7
                props=Tuple((AC.Energy(),AC.Forces(),AC.Virial())[j] for j=1:3 if (subset>>(j-1))&1==1)
                result=GA.gpu_result(GA.execute_gpu!(ex;properties=props))
                @test length(result)==length(props)
                for key in keys(result)
                    actual=key===:energy ? only(result[key]) : result[key]
                    @test actual ≈ ref[key] atol=tol rtol=T===Float64 ? 1e-10 : 5e-4
                    @test eltype(result[key])===T
                end
            end
            owned=GA.execute_gpu!(ex)
            old=GA.gpu_result(owned)
            old_density=ex.U; old_connectivity=ex.workspace.cache.neighbors
            GA.execute_gpu!(ex;owning=false)
            @test GA.gpu_result(owned).forces==old.forces
            @test ex.U===old_density && ex.workspace.cache.neighbors===old_connectivity
            @test old.forces ≈ ref.forces atol=tol rtol=T===Float64 ? 1e-10 : 5e-4
        end
        for centers in (Int[],[3,1])
            ex=GA.GPUExecutor(sys,b.model,b.parameters;centers,batch_size=1)
            y=GA.gpu_result(GA.execute_gpu!(ex)); r=cpu_reference(sys,cpu.model,cpu.parameters;centers)
            @test only(y.energy) ≈ r.energy atol=tol
            @test y.forces ≈ r.forces atol=tol rtol=5e-4
            @test y.virial ≈ r.virial atol=tol rtol=5e-4
        end
        for n in (0,1)
            tiny=TestSystem(sys.xyz[:,1:n],sys.H,sys.pbc,sys.labels[1:n]); ex=GA.GPUExecutor(tiny,b.model,b.parameters)
            y=GA.gpu_result(GA.execute_gpu!(ex)); r=cpu_reference(tiny,cpu.model,cpu.parameters)
            @test only(y.energy) ≈ r.energy atol=tol
            @test y.forces==r.forces
            @test y.virial==r.virial
        end
    end
end
@testset "GPU ACE derivatives, updates, and custom radials" begin
    m,p,_=GA.load_resolved(joinpath(@__DIR__,"../data/ace_ref/resolved_model_v1.toml"))
    sys=TestSystem([.4 2.6 1.1;0. .2 2.8;0. .3 .1],Matrix{Float64}(I,3,3)*16,(false,false,false),[:Ti,:Al,:Ti])
    b=gpu_bundle(m,p,Float64)
    ex=GA.GPUExecutor(sys,b.model,b.parameters;batch_size=2)
    y=GA.gpu_result(GA.execute_gpu!(ex)); saved=deepcopy(y)
    function evaluate_system(s)
        GA.prepare_executor!(ex,s)
        GA.gpu_result(GA.execute_gpu!(ex))
    end
    h=1e-5
    for a=1:3,i=1:3
        plus=deepcopy(sys); minus=deepcopy(sys); plus.xyz[a,i]+=h; minus.xyz[a,i]-=h
        fd=-(only(evaluate_system(plus).energy)-only(evaluate_system(minus).energy))/(2h)
        @test y.forces[a,i] ≈ fd atol=2e-7 rtol=2e-7
    end
    for a=1:3,c=1:3
        strain=zeros(3,3); strain[c,a]=h
        plus=TestSystem((I+strain)*sys.xyz,sys.H,sys.pbc,sys.labels)
        minus=TestSystem((I-strain)*sys.xyz,sys.H,sys.pbc,sys.labels)
        fd=-(only(evaluate_system(plus).energy)-only(evaluate_system(minus).energy))/(2h)
        @test y.virial[a,c] ≈ fd atol=2e-7 rtol=2e-7
    end
    GA.prepare_executor!(ex,sys)
    changed=deepcopy(p); changed.WB .*= 1.1
    converted=GA.adapt_model(m,changed;array_type=array_constructor)
    updated=GA.gpu_result(GA.execute_gpu!(ex;parameters=converted.parameters))
    @test updated.forces ≈ cpu_reference(sys,m,changed).forces atol=1e-9
    @test y.forces==saved.forces
    for n in (5,1,0,3)
        xyz=hcat(sys.xyz,[8. 10.;0. .1;0. 0.])[:,1:n]
        s=TestSystem(xyz,sys.H,sys.pbc,[:Ti,:Al,:Ti,:Al,:Ti][1:n])
        GA.prepare_executor!(ex,s)
        @test GA.gpu_result(GA.execute_gpu!(ex)).forces ≈ cpu_reference(s,m,changed).forces atol=1e-9
    end
    # An explicit custom analytic radial shares the public callable protocol.
    meta=GA.RadialMetadata(m.radial_spec;cutoff=5.5,origin=:regular)
    channel=GA.CallableRadial(GA.SmoothExponential(),meta;config=(cutoff=5.5,powers=collect(0:length(m.radial_spec)-1)))
    replacement=GA.replace_radial(m,p,fill(channel,2,2);parameters=(rate=[.3],amplitude=[1.2]))
    custom=gpu_bundle(replacement.model,replacement.parameters,Float64)
    cx=GA.GPUExecutor(sys,custom.model,custom.parameters;strategy=:compact)
    cr=GA.gpu_result(GA.execute_gpu!(cx)); ref=cpu_reference(sys,replacement.model,replacement.parameters)
    @test only(cr.energy) ≈ ref.energy atol=1e-9
    @test cr.forces ≈ ref.forces atol=1e-9
    @test cr.virial ≈ ref.virial atol=1e-9
    scalar_channel=GA.CallableRadial(ScalarParameterRadial(),meta;config=(cutoff=5.5,))
    scalar=GA.replace_radial(m,p,fill(scalar_channel,2,2);parameters=(rate=.3,))
    sb=gpu_bundle(scalar.model,scalar.parameters,Float64)
    sx=GA.GPUExecutor(sys,sb.model,sb.parameters)
    previous=GA.gpu_result(GA.execute_gpu!(sx))
    nextps=merge(sx.parameters,(custom_radial=(rate=.6,),))
    changed_scalar=GA.gpu_result(GA.execute_gpu!(sx;parameters=nextps))
    reference_scalar=cpu_reference(sys,scalar.model,merge(scalar.parameters,(custom_radial=(rate=.6,),)))
    @test sx.parameters.custom_radial.rate==.6
    @test changed_scalar.forces ≈ reference_scalar.forces atol=1e-9
    @test maximum(abs,changed_scalar.forces-previous.forces)>1e-10
    coincident=deepcopy(sys); coincident.xyz[:,2]=coincident.xyz[:,1]
    GA.prepare_executor!(ex,coincident)
    @test_throws DomainError GA.execute_gpu!(ex;properties=(AC.Forces(),))
end
@testset "GPU ACE geometry, filters, and invariance" begin
    m,p,_=GA.load_resolved(joinpath(@__DIR__,"../data/ace_ref/resolved_model_v1.toml"))
    sys=TestSystem([15.8 2.1 1.1;0. .2 2.8;0. .3 .1],Matrix{Float64}(I,3,3)*18,(true,false,true),[:Ti,:Al,:Ti])
    for model in (m,GA.project_species(m,[:Al],0:1),GA.mask_basis(m,GA.BasisSelection(m;indices=[1,5,20])))
        b=gpu_bundle(model,p,Float64)
        ex=GA.GPUExecutor(sys,b.model,b.parameters;config=GN.GPUConfig(skin=.6))
        for step in (0.,.02,.12,.35)
            moved=deepcopy(sys); moved.xyz[1,1]+=step
            copyto!(ex.positions,array_constructor(moved.xyz))
            result=GA.gpu_result(GA.execute_gpu!(ex))
            ref=cpu_reference(moved,model,p)
            @test only(result.energy) ≈ ref.energy atol=1e-9
            @test result.forces ≈ ref.forces atol=1e-9
            @test result.virial ≈ ref.virial atol=1e-9
        end
        @test ex.workspace.rebuilds>=2
    end
    base=TestSystem([.4 2.6 1.1;0. .2 2.8;0. .3 .1],sys.H,(false,false,false),sys.labels)
    b=gpu_bundle(m,p,Float64)
    evaluate_case(s)=GA.gpu_result(GA.execute_gpu!(GA.GPUExecutor(s,b.model,b.parameters)))
    ref=evaluate_case(base)
    translated=TestSystem(base.xyz .+ [3.,-2.,1.],base.H,base.pbc,base.labels)
    @test evaluate_case(translated).forces ≈ ref.forces atol=1e-9
    perm=[3,1,2]; permuted=TestSystem(base.xyz[:,perm],base.H,base.pbc,base.labels[perm])
    @test evaluate_case(permuted).forces ≈ ref.forces[:,perm] atol=1e-9
    rot=[cos(.3) -sin(.3) 0.;sin(.3) cos(.3) 0.;0. 0. 1.]
    rotated=evaluate_case(TestSystem(rot*base.xyz,base.H,base.pbc,base.labels))
    @test only(rotated.energy) ≈ only(ref.energy) atol=1e-9
    @test rotated.forces ≈ rot*ref.forces atol=1e-9
    @test rotated.virial ≈ rot*ref.virial*rot' atol=1e-9
    @test vec(sum(ref.forces;dims=2)) ≈ zeros(3) atol=1e-12
end
if "--amdgpu" in ARGS
@testset "AMD selected-device placement and calculator boundary" begin
    @test !isnothing(Base.get_extension(PotentialFrameWork,:PotentialFrameWorkAMDGPUExt))
    m,p,_=GA.load_resolved(joinpath(@__DIR__,"../data/ace_ref/resolved_model_v1.toml"))
    sys=TestSystem([.4 2.6 1.1;0. .2 2.8;0. .3 .1],Matrix{Float64}(I,3,3)*16,(false,false,false),[:Ti,:Al,:Ti])
    b=gpu_bundle(m,p,Float32)
    @test_throws DimensionMismatch GA.GPUExecutor(sys,b.model,(WB=AMDGPU.ROCArray(zeros(Float32,1,2)),Wpair=b.parameters.Wpair))
    b64=gpu_bundle(m,p,Float64)
    @test GA.automatic_gpu_strategy(b64.model,b64.model.tables.reference,127)===:staged
    @test GA.automatic_gpu_strategy(b64.model,b64.model.tables.reference,128)===:compact
    @test GA.automatic_gpu_strategy(b64.model,b64.model.tables.reference,1024)===:compact
    @test GA.automatic_gpu_strategy(b64.model,b64.model.tables.reference,1025)===:staged
    @test GA.automatic_gpu_strategy(b.model,b.model.tables.reference,128)===:staged
    resident=GA.GPUExecutor(sys,b.model,b.parameters)
    current=AMDGPU.device()
    other=filter(!=(current),AMDGPU.devices())
    if !isempty(other)
        try
            AMDGPU.device!(first(other))
            @test_throws ArgumentError GA.GPUExecutor(sys,b.model,b.parameters)
            @test_throws ArgumentError GA.prepare_executor!(resident,sys)
        finally
            AMDGPU.device!(current)
        end
    end
    calc=GA.ACECalculator(m,p;array_type=AMDGPU.ROCArray,precision=Float32)
    y=AC.energy_forces_virial(sys,calc)
    @test length(y.forces)==3
    mktemp() do path,io
        close(io); GA.save_calculator(path,calc)
        loaded=GA.ACECalculator(path;array_type=AMDGPU.ROCArray,precision=Float32)
        @test reduce(hcat,[ustrip.(u"eV/Å",v) for v in AC.energy_forces_virial(sys,loaded).forces]) ≈ reduce(hcat,[ustrip.(u"eV/Å",v) for v in y.forces])
        @test loaded.model.source.coefficient_precision===Float32
    end
end
end
@testset "GPU ACE constant angular and empty invariant bases" begin
    m,p,_=GA.load_resolved(joinpath(@__DIR__,"../data/ace_ref/resolved_model_v1.toml"))
    constant=GA.ResolvedModel(copy(m.atomic_numbers),zeros(2),m.radial,m.pair_radial,m.radial_spec,m.pair_spec,0,
        [(1,1)],[[1],[1,1]],GA.sparse([1,2],[1,2],[1.,1.],2,2),m.source_sha256)
    cp=(WB=[.3 -.2;.7 .4],Wpair=p.Wpair)
    empty=GA.prune_basis(m,p,GA.BasisSelection(m;indices=Int[]))
    sys=TestSystem([.4 2.6 1.1;0. .2 2.8;0. .3 .1],Matrix{Float64}(I,3,3)*16,(true,true,true),[:Ti,:Al,:Ti])
    for (model,ps) in ((constant,cp),(empty.model,empty.parameters)),T in (Float64,Float32)
        b=gpu_bundle(model,ps,T); cpu=GA.adapt_model(model,ps;precision=T)
        ex=GA.GPUExecutor(sys,b.model,b.parameters;batch_size=1)
        result=GA.gpu_result(GA.execute_gpu!(ex)); ref=cpu_reference(sys,cpu.model,cpu.parameters)
        @test only(result.energy) ≈ ref.energy rtol=1e-5
        @test result.forces ≈ ref.forces atol=T===Float32 ? 1e-4 : 1e-10
        @test result.virial ≈ ref.virial atol=T===Float32 ? 1e-4 : 1e-10
        @test isempty(ex.V) && isempty(ex.dV)
    end
end
