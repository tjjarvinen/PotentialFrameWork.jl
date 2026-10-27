# Requires AMDGPU, AtomsBuilder, AtomsSystems in an application environment.
using PotentialFrameWork, AMDGPU, AtomsBase, AtomsBuilder, AtomsSystems
using Unitful, Random, Statistics, TOML, Pkg
import AtomsCalculators as AC
const NL=PotentialFrameWork.NeighbourList
const PP=PotentialFrameWork.PairPotentials
AMDGPU.functional() || error("A real AMDGPU device is required")
AMDGPU.device!(only(filter(d->occursin("6800 XT",string(d)),AMDGPU.devices())))
const backend=NL.KA.get_backend(ROCArray(zeros(1)))
const repetitions=parse(Int,get(ENV,"PFW_GPU_REPETITIONS","7"))
const rows=Dict{String,Any}[]
function measure(f)
    f(); NL.KA.synchronize(backend)
    times=Float64[]
    for _=1:repetitions
        NL.KA.synchronize(backend)
        t=time_ns(); f(); NL.KA.synchronize(backend)
        push!(times,(time_ns()-t)*1e-9)
    end
    Dict("median"=>median(times),"minimum"=>minimum(times),"maximum"=>maximum(times))
end
function structure(m,kind)
    Random.seed!(20261027+m)
    built=repeat(AtomsBuilder.bulk(:Ar;cubic=true,a=5.26u"Å"),(m,m,m))
    AtomsBuilder.rattle!(built,.03u"Å")
    sys=AtomsSystems.generic_system(built)
    if kind!=:argon
        labels=[isodd(i) ? :Ar : :Ne for i=1:length(sys)]
        sys=AtomsSystems.generic_system(labels,collect(position(sys,:));
            cell_vectors=cell_vectors(sys),periodicity=kind===:slab ? (true,true,false) : (true,true,true))
    end
    sys
end
function stage_profile(input,request)
    w=input.workspace; n=length(input.identities); B=w.config.workgroup
    stages=Dict{String,Any}()
    stages["key_generation"]=Dict(pairs(measure(()->NL.gpu_records!(backend,B)(w.records,input.positions,input.identities,w.lengths,w.periodic,w.counts,w.radius;ndrange=n))))
    stages["keys_and_sort"]=Dict(pairs(measure(()->begin
        NL.gpu_records!(backend,B)(w.records,input.positions,input.identities,w.lengths,w.periodic,w.counts,w.radius;ndrange=n)
        NL.gpu_sort!(w.records,w.sorttemp)
    end))) # includes key generation to sort identical unsorted records each sample
    stages["unpack_ranges"]=Dict(pairs(measure(()->NL.gpu_unpack!(backend,B)(w.x,w.ends,w.permutation,w.species,w.records,input.labels,n;ndrange=n))))
    args=(w.force,w.partials,w.x,w.records,w.ends,w.species,w.permutation,w.lengths,w.periodic,w.counts,w.radius,input.prepared,input.parameters,PP.properties(request),n)
    run= w.config.strategy===:central ? ()->NL.gpu_central!(backend,B)(args...;ndrange=n) :
        ()->NL.gpu_candidate!(backend,B)(args...,Val(B),Val(w.config.tile);ndrange=cld(n,w.config.tile)*B)
    stages["interaction"]=Dict(pairs(measure(run)))
    stages["reduction"]=Dict(pairs(measure(()->NL.gpu_reduce!(w.totals,w.partials,w.reducetemp))))
    stages
end
scratch_bytes(x::AbstractArray)=sizeof(x)
scratch_bytes(x::NamedTuple)=sum(scratch_bytes,values(x);init=0)
scratch_bytes(::Nothing)=0
scratch_bytes(x)=scratch_bytes(x.buffers)
for m in (14,29),kind in (:argon,:mixed,:slab)
    sys=structure(m,kind)
    p=kind===:argon ? PP.lennard_jones() : SmoothPair([:Ar,:Ne],[.0103 .008;.008 .006],8.5)
    original_reference=evaluate(sys,p;precision=Float64,width=4,parallel=true)
    for T in (Float32,Float64), strategy in (:central,:candidate), request in (EnergyOnly(),EnergyForcesVirial())
        NL.KA.synchronize(backend); prep_start=time_ns()
        input=PP.prepare_gpu_execution(sys,p,backend;precision=T,config=NL.GPUConfig(;strategy))
        NL.KA.synchronize(backend); initial_prepare_seconds=(time_ns()-prep_start)*1e-9
        first_start=time_ns(); PP.execute_gpu!(input;request); NL.KA.synchronize(backend)
        first_execution_seconds=(time_ns()-first_start)*1e-9 # compilation/capacity growth if needed
        # Numerical correctness is a benchmark gate, including original-order forces.
        y=PP.evaluate_gpu!(input.workspace,sys,p;request)
        ref=evaluate(sys,p;precision=T,width=4,request,parallel=true)
        tol=T===Float32 ? 2e-4 : 2e-10
        isapprox(ustrip(y.energy),ustrip(ref.energy);rtol=tol,atol=tol) || error("energy mismatch")
        ferr=haskey(y,:forces) ? maximum(abs,ustrip.(y.forces-ref.forces)) : 0.0
        ferr<tol || error("force mismatch: $ferr")
        if haskey(y,:virial)
            isapprox(ustrip.(y.virial),ustrip.(ref.virial);atol=tol,rtol=tol) || error("virial mismatch")
        end
        input=PP.prepare_gpu_execution(input.workspace,sys,p)
        numeric=measure(()->PP.execute_gpu!(input;request))
        public=measure(()->PP.evaluate_gpu!(input.workspace,sys,p;request))
        cpu=measure(()->evaluate(sys,p;precision=T,width=4,request,parallel=true))
        profile=stage_profile(input,request)
        profile["extraction"]=measure(()->PP.gpu_inputs(sys,p,T,input.workspace.config))
        host=PP.gpu_inputs(sys,p,T,input.workspace.config)
        profile["geometry_upload"]=measure(()->begin
            copyto!(input.positions,host.positions); copyto!(input.identities,host.identities)
        end)
        profile["host_output_copy"]=measure(()->begin
            Array(input.workspace.totals)
            request isa EnergyOnly || Array(input.workspace.force)
        end)
        input=PP.prepare_gpu_execution(input.workspace,sys,p)
        result=PP.execute_gpu!(input;request); NL.KA.synchronize(backend)
        host_allocated=@allocated begin
            PP.execute_gpu!(input;request)
            NL.KA.synchronize(backend)
        end
        row=Dict{String,Any}("n"=>length(sys),"replication"=>m,"kind"=>string(kind),"precision"=>string(T),
            "strategy"=>string(strategy),"request"=>string(typeof(request)),"skin"=>0.,"seed"=>20261027+m,
            "perturbation_A"=>.03,"cutoff_A"=>8.5,"max_force_error_eV_A"=>ferr,
            "max_atom_force_error_eV_A"=>(haskey(y,:forces) ? maximum(sqrt.(sum(abs2,ustrip.(y.forces-ref.forces);dims=1))) : 0.),
            "max_atom_force_error_original_eV_A"=>(haskey(y,:forces) ? maximum(sqrt.(sum(abs2,ustrip.(y.forces)-ustrip.(original_reference.forces);dims=1))) : 0.),
            "max_force_error_original_eV_A"=>(haskey(y,:forces) ? maximum(abs,ustrip.(y.forces)-ustrip.(original_reference.forces)) : 0.),
            "energy_error_original_per_atom_eV"=>abs(ustrip(y.energy)-ustrip(original_reference.energy))/length(sys),
            "energy_error_per_atom_eV"=>abs(ustrip(y.energy-ref.energy))/length(sys),
            "numeric"=>Dict(pairs(numeric)),"numeric_atoms_per_second"=>length(sys)/numeric["median"],"public"=>Dict(pairs(public)),"cpu"=>Dict(pairs(cpu)),
            "main_buffer_bytes"=>sum(sizeof,(input.workspace.records,input.workspace.x,input.workspace.ends,input.workspace.species,input.workspace.permutation,input.workspace.partials,input.workspace.totals,input.workspace.force)),
            "upload_buffer_bytes"=>sum(a->isnothing(a) ? 0 : sizeof(a),input.workspace.transfer.arrays),
            "scratch_bytes"=>scratch_bytes(input.workspace.sorttemp)+scratch_bytes(input.workspace.reducetemp),
            "owning_result_bytes"=>sum(sizeof,values(result)),"host_allocated_bytes"=>host_allocated,
            "cell_lengths_A"=>[ustrip(u"Å",cell_vectors(sys)[a][a]) for a=1:3],
            "periodicity"=>collect(periodicity(sys)),"density_atoms_A3"=>length(sys)/(5.26*m)^3,
            "species_fractions"=>(kind===:argon ? [1.] : [.5,.5]),
            "workgroup"=>64,"central_tile"=>4,"initial_prepare_seconds"=>initial_prepare_seconds,"first_execution_seconds"=>first_execution_seconds,"stages"=>profile)
        push!(rows,row); println(row)
        open(joinpath(@__DIR__,"gpu_pair_measurements.toml"),"w") do io
            TOML.print(io,Dict("julia"=>string(VERSION),"cpu"=>Sys.CPU_NAME,"device"=>string(AMDGPU.device()),
                "repetitions"=>repetitions,"julia_threads"=>Threads.nthreads(),"dependencies"=>Dict(d.name=>string(d.version) for d in values(Pkg.dependencies()) if d.version!==nothing),"measurements"=>rows))
        end
    end
end
# Evolving trajectory: retain raw coordinates, accumulate small bounded displacements.
for strategy in (:central,:candidate),skin in (0.,.5)
    sys=structure(14,:argon); p=PP.lennard_jones()
    input=PP.prepare_gpu_execution(sys,p,backend;precision=Float32,config=NL.GPUConfig(;strategy,skin))
    PP.execute_gpu!(input); PP.execute_gpu!(input); NL.KA.synchronize(backend)
    x=Array(input.positions); reference=copy(x); timings=Float64[]; rebuilt=Bool[]
    velocity=[Float32(.005*sin(a+0.73*i)) for a=1:3,i=1:length(sys)]
    for step=1:40
        x.=reference.+Float32(step).*velocity
        copyto!(input.positions,x)
        previous_rebuilds=input.workspace.rebuilds
        NL.KA.synchronize(backend); start=time_ns()
        result=PP.execute_gpu!(input); NL.KA.synchronize(backend)
        push!(timings,(time_ns()-start)*1e-9); push!(rebuilt,input.workspace.rebuilds>previous_rebuilds)
        # Compare evolving coordinates to fresh GPU, outside timed region.
        fresh=PP.prepare_gpu_execution(sys,p,backend;precision=Float32,config=NL.GPUConfig(;strategy))
        copyto!(fresh.positions,x); ref=PP.execute_gpu!(fresh); NL.KA.synchronize(backend)
        isapprox(Array(result.forces),Array(ref.forces);atol=3e-5,rtol=3e-5) || error("skin trajectory mismatch")
    end
    row=Dict("strategy"=>string(strategy),"skin"=>skin,"steps"=>40,"rebuild_count"=>count(rebuilt),
        "precision"=>"Float32","n"=>length(sys),"step_amplitude_A"=>.005,
        "skin_connectivity_bytes"=>(skin==0 ? 0 : sum(sizeof,(input.workspace.cache.offsets,input.workspace.cache.neighbors,input.workspace.cache.reference,input.workspace.cache.distance2,input.workspace.cache.identities,input.workspace.cache.labels))+scratch_bytes(input.workspace.cache.scanworkspace)+scratch_bytes(input.workspace.cache.maxworkspace)),
        "amortized_seconds"=>mean(timings),"rebuild_seconds"=>timings[rebuilt],"reuse_seconds"=>timings[.!rebuilt])
    open(joinpath(@__DIR__,"gpu_skin_$(strategy)_$(skin).toml"),"w") do io; TOML.print(io,row); end
end
