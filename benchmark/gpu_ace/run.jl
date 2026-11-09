# Explicit benchmark environment; see README.md for setup and timed boundaries.
using PotentialFrameWork, BenchmarkTools, LinearAlgebra, TOML, SHA, Statistics
const A=PotentialFrameWork.ACE
const NL=PotentialFrameWork.NeighbourList
const AB=A.AtomsBase
const U=A.Unitful
const AC=A.AC
function options(args)
    d=Dict("preset"=>"quick","backend"=>"cpu","models"=>"tial,ti,tialv","precision"=>"Float64","strategy"=>"staged",
           "samples"=>"5","seconds"=>"0.3","output"=>"gpu-ace-report.toml","data"=>"","sizes"=>"","parallel"=>"true","spacing"=>"3.4")
    iseven(length(args)) || error("Use --option value pairs")
    for k=1:2:length(args)
        key=replace(args[k],r"^--"=>""); haskey(d,key) || error("unknown option $(args[k])"); d[key]=args[k+1]
    end
    d["preset"] in ("quick","full") || error("preset must be quick or full")
    d
end
include("workloads.jl")
include("provenance.jl")
function structure(g,z;spacing=3.4)
    data=bcc_geometry(g,z;spacing); X=data.positions; H=data.cell; n=size(X,2)
    cells=ntuple(a->A.SVector{3}(H[:,a])*U.u"Å",3)
    AB.atomic_system([AB.ChemicalSpecies(data.species[i])=>A.SVector{3}(X[:,i])*U.u"Å" for i=1:n],cells,(true,true,true))
end

function stats(trial,n)
    times=trial.times*1e-9
    Dict("median_s"=>median(times),"minimum_s"=>minimum(times),"maximum_s"=>maximum(times),"std_s"=>length(times)>1 ? std(times) : 0.,
         "samples"=>length(times),"raw_seconds"=>times,"atoms_per_second"=>n/median(times),"host_bytes"=>trial.memory,"host_allocations"=>trial.allocs,"raw_gc_ns"=>trial.gctimes,"device_allocations"=>"unavailable")
end
function measure(f,n,samples,seconds,sync;reset=()->nothing)
    reset(); f(); sync()
    completed()=begin f(); sync(); nothing end
    benchmark=@benchmarkable $completed() setup=($reset()) evals=1
    trial=run(benchmark;samples,seconds)
    stats(trial,n)
end
function accuracy(y,ref)
    ev=abs((y.energy isa Number ? y.energy : only(y.energy))-ref.energy)
    fd=Float64.(y.forces)-Float64.(ref.forces); wd=Float64.(y.virial)-Float64.(ref.virial)
    ef=maximum(abs,fd;init=0.); ew=maximum(abs,wd;init=0.)
    all(isfinite,(ev,ef,ew)) || error("nonfinite numerical check")
    Dict("absolute_energy_eV"=>ev,"energy_per_atom_eV"=>size(ref.forces,2)==0 ? 0. : ev/size(ref.forces,2),
         "force_max_component_eV_per_A"=>ef,"virial_max_component_eV"=>ew,
         "force_rms_component_eV_per_A"=>isempty(fd) ? 0. : sqrt(sum(abs2,fd)/length(fd)),
         "force_max_atom_vector_eV_per_A"=>maximum(sqrt.(vec(sum(abs2,fd;dims=1)));init=0.),
         "virial_rms_component_eV"=>sqrt(sum(abs2,wd)/9))
end
function trajectory(sys,model,ps,storage,T,strategy,n,samples,seconds)
    inp=A.ace_inputs(sys,model)
    # Replay exactly the same deterministic sequence; no continuously drifting input.
    shifts=(0.,.025,.075,.16,.4)
    Xs=[T.(inp.positions .+ shift.*[sin(i*(a+.31)) for a=1:3,i=1:n]) for shift in shifts]
    if storage===Array
        bundle=A.adapt_model(model,ps;precision=T)
        env=A.FullEnvironment(width=4)
        frames=[AB.atomic_system([AB.species(sys,i)=>A.SVector{3}(X[:,i])*U.u"Å" for i=1:n],AB.cell_vectors(sys),AB.periodicity(sys)) for X in Xs]
        ex=A.CPUExecutor(bundle.model,env;precision=T,angular_backend=:block)
        replay=()->begin
            for frame in frames
                input=A.ace_inputs(A.precision_geometry(frame,T),bundle.model;environment=env,precision=T)
                A.prepare_executor!(ex,env)
                A.execute_prepared!(ex,T.(input.positions),T.(input.cell),bundle.parameters;parallel=true)
            end
        end
        result=measure(replay,n*length(shifts),samples,seconds,()->nothing)
        result["boundary"]="Five complete CPU extraction/topology/prepared steps; same sequence each sample"
    else
        bundle=A.adapt_model(model,ps;array_type=storage,precision=T)
        skin=min(.6,.5*(minimum(inp.cell[a,a] for a=1:3)/2-A.NI.search_cutoff(model)))
        ex=A.GPUExecutor(sys,bundle.model,bundle.parameters;strategy,config=NL.GPUConfig(;skin))
        device_frames=map(storage,Xs)
        replay=()->begin
            for X in device_frames
                copyto!(ex.positions,X); A.execute_gpu!(ex;owning=false)
            end
        end
        reset=()->begin ex.workspace.cache=nothing; ex.workspace.rebuilds=0; nothing end
        result=measure(replay,n*length(shifts),samples,seconds,()->NL.KA.synchronize(ex.workspace.backend);reset)
        result["rebuilds_last_sample"]=ex.workspace.rebuilds
        result["skin_A"]=skin
        result["boundary"]="Five resident device-to-device position copies and evaluations; cache cleared in sample setup; first rebuild included"
        result["sequence_device_bytes"]=sum(sizeof,device_frames)
    end
    result["steps"]=length(shifts); result
end

function main()
    start=time(); opts=options(ARGS); BLAS.set_num_threads(1)
    storage=if opts["backend"]=="amd"
        AMDGPU.functional() || error("AMD requested but unavailable")
        AMDGPU.device!(only(filter(d->occursin("6800 XT",string(d)),AMDGPU.devices())))
        AMDGPU.ROCArray
    elseif opts["backend"]=="cpu"
        Array
    else
        error("backend must be cpu or amd")
    end
    T=opts["precision"]=="Float32" ? Float32 : opts["precision"]=="Float64" ? Float64 : error("precision")
    samples=parse(Int,opts["samples"]); seconds=parse(Float64,opts["seconds"])
    sizes=isempty(opts["sizes"]) ? (opts["preset"]=="quick" ? [4,18] : [4,8,18,24]) : parse.(Int,split(opts["sizes"],','))
    records=Dict[]; omitted=String[]
    for name in split(opts["models"],',')
        path=name=="tial" ? joinpath(@__DIR__,"../../data/ace_ref/resolved_model_v1.toml") : joinpath(opts["data"],name*"-resolved.toml")
        if !isfile(path); push!(omitted,"$name: no independently verified resolved export at $path"); continue; end
        if name!="tial"
            verification=joinpath(opts["data"],"verification-cpu.toml")
            isfile(verification) || error("Run test/gpu_ace_public.jl on CPU before benchmarking public exports")
            TOML.parsefile(verification)["resolved_sha256"][name]==bytes2hex(sha256(read(path))) || error("public export changed since independent verification")
        end
        model,ps,_=A.load_resolved(path); bundle=A.adapt_model(model,ps;array_type=storage,precision=T)
        for g in (opts["preset"]=="quick" && name!="tial" ? sizes[1:1] : sizes)
            g>=4 || error("periodic single-image workloads need replication >=4")
            sys=structure(g,model.atomic_numbers;spacing=parse(Float64,opts["spacing"])); n=length(sys)
            # Common Float64 reference uses existing CPU mathematical executor.
            ref=A.precision_evaluate(sys,(;model,parameters=ps),(AC.Energy(),AC.Forces(),AC.Virial()))
            source_geometry=A.precision_geometry(sys,Float64)
            input=A.ace_inputs(storage===Array ? A.precision_geometry(sys,T) : sys,bundle.model isa A.GPUACEModel ? bundle.model.source : bundle.model;precision=T)
            if storage===Array
                ex=A.CPUExecutor(bundle.model,input.environment;precision=T,angular_backend=:block)
                X=T.(input.positions); H=T.(input.cell)
                prepared=()->A.execute_prepared!(ex,X,H,bundle.parameters;parallel=opts["parallel"]=="true")
                sync=()->nothing
                prepared(); actual=A.prepared_result(ex)
                memory=Dict{String,Any}("site_scratch_bytes"=>Base.summarysize(ex.scratch),"outputs_bytes"=>sizeof(ex.force)+sizeof(ex.virials))
            else
                ex=A.GPUExecutor(sys,bundle.model,bundle.parameters;strategy=Symbol(opts["strategy"]))
                prepared=()->A.execute_gpu!(ex;owning=false)
                sync=()->NL.KA.synchronize(ex.workspace.backend)
                actual=A.gpu_result(prepared())
                stages=Dict{Symbol,Float64}(); A.execute_gpu!(ex;owning=false,diagnostics=stages)
                memory=Dict{String,Any}("density_bytes"=>sizeof(ex.U)+sizeof(ex.V),"adjoint_bytes"=>sizeof(ex.dU)+sizeof(ex.dV),
                    "weights_bytes"=>sizeof(ex.weights),"output_force_bytes"=>sizeof(ex.workspace.force),
                    "connectivity_bytes"=>sizeof(ex.workspace.cache.neighbors)+sizeof(ex.workspace.cache.offsets),
                    "model_host_descriptor_bytes"=>Base.summarysize(model),"model_device_array_bytes"=>A.gpu_memory(ex).model_tables)
            end
            if storage!==Array
                memory["density_shared_payload_bytes_per_group"]=ex.strategy===:baseline ? 0 : 2*8*8*sizeof(T)
                memory["compact_shared_payload_bytes_per_group"]=ex.strategy===:compact ? 13*4*16*sizeof(T) : 0
                memory["compiler_register_and_spill_bytes"]="unavailable"
                memory["resident_array_payloads"]=Dict(string(k)=>v for (k,v) in pairs(A.gpu_memory(ex)))
            end
            directed_edges=storage===Array ? count(input.environment.active) : only(Array(view(ex.workspace.cache.offsets,n:n)))
            stage_report=storage===Array ? Dict{String,Float64}() : Dict(string(k)=>v for (k,v) in stages)
            errors=accuracy(actual,ref)
            oracle_timing=Dict()
            oracle_path=joinpath(opts["data"],name*"-oracle.toml")
            if isfile(oracle_path)
                oracle=TOML.parsefile(oracle_path)
                oracle["source_sha256"]==model.source_sha256 || error("oracle model hash mismatch")
                for r in oracle["records"]
                    if length(r["species"])==n && maximum(abs,reshape(r["positions"],3,n)-source_geometry.positions)<=1e-12
                        isapprox(ref.energy,r["energy"];atol=2e-7,rtol=1e-10) || error("independent oracle energy mismatch")
                        isapprox(ref.virial,reshape(r["virial"],3,3);atol=2e-6,rtol=1e-9) || error("independent oracle virial mismatch")
                        isapprox(ref.forces,reshape(r["forces"],3,n);atol=2e-7,rtol=1e-9) || error("independent oracle force mismatch")
                        oracle_timing=Dict("median_s"=>get(r,"median_s",NaN),"raw_ns"=>get(r,"raw_ns",Float64[]),"boundary"=>get(r,"timed_boundary","not timed"),
                            "julia"=>oracle["julia"],"ACEpotentials"=>oracle["ACEpotentials"],"precision"=>"Float64 original")
                    end
                end
            end
            # Implementation precision checks, separate from application sensitivity.
            if T===Float64
                errors["energy_per_atom_eV"]<=1e-9 || error("Float64 energy check failed")
                errors["force_max_component_eV_per_A"]<=1e-7 || error("Float64 force check failed")
                errors["virial_max_component_eV"]<=1e-6 || error("Float64 virial check failed")
            end
            if T===Float32
                rounded=A.precision_evaluate(sys,A.adapt_model(model,ps;precision=Float32),(AC.Energy(),AC.Forces(),AC.Virial()))
                actual_energy=actual.energy isa Number ? actual.energy : only(actual.energy)
                abs(actual_energy-rounded.energy)<=16eps(Float32)*max(abs(rounded.energy),n) || error("Float32 energy implementation check failed")
                errors["implementation_energy_abs_eV"]=abs(actual_energy-rounded.energy)
                errors["implementation_virial_max_component_eV"]=maximum(abs,actual.virial-rounded.virial;init=0.)
                isapprox(actual.forces,rounded.forces;atol=2e-3,rtol=1e-3) || error("Float32 implementation check failed")
                isapprox(actual.virial,rounded.virial;atol=2e-2,rtol=1e-3) || error("Float32 virial implementation check failed")
                errors["implementation_force_max_component_eV_per_A"]=maximum(abs,actual.forces-rounded.forces;init=0.)
            end
            result=measure(prepared,n,samples,seconds,sync)
            calculator=A.ACECalculator(model,ps;array_type=storage,precision=T,angular_backend=:block)
            storage===Array || AC.energy_forces_virial(sys,calculator;strategy=Symbol(opts["strategy"]))
            complete_call=storage===Array ? (()->AC.energy_forces_virial(sys,calculator;parallel=opts["parallel"]=="true")) : (()->AC.energy_forces_virial(sys,calculator))
            complete=measure(complete_call,n,samples,seconds,sync)
            sequence=opts["preset"]=="full" && g==first(sizes) ? trajectory(sys,model,ps,storage,T,Symbol(opts["strategy"]),n,samples,seconds) : Dict()
            push!(records,Dict("model"=>name,"model_sha256"=>model.source_sha256,"resolved_sha256"=>bytes2hex(sha256(read(path))),
                "atoms"=>n,"replication"=>g,"structure_sha256"=>bytes2hex(sha256(string(source_geometry.positions,source_geometry.cell,model.atomic_numbers))),"position_rounding_max_A"=>maximum(abs,Float64.(T.(source_geometry.positions))-source_geometry.positions;init=0.),"precision"=>string(T),"requested_strategy"=>opts["strategy"],"strategy"=>storage===Array ? "colored CPU" : string(ex.strategy),"density_components"=>length(model.density),"species_atomic_numbers"=>model.atomic_numbers,"species_counts"=>[count(i->Int(AB.atomic_number(AB.species(sys,i)))==z,1:n) for z in model.atomic_numbers],
                "directed_edges"=>directed_edges,"complete_strategy"=>storage===Array ? "colored CPU" : string(calculator.executor.strategy),
                "products"=>length(model.products),"readout_coordinates"=>size(model.coupling,1),"readout_coordinate_kind"=>string(model.coordinate_kind),
                "radial_channels"=>length(model.radial_spec),"angular_channels"=>(model.lmax+1)^2,
                "published_invariant_basis_count"=>model.source_sha256=="aa4c9f27bbb17709ce31ab0bd7786b9f363b458f6e0ed2cfcbfe9a7a57834e8b" ? 1809 : model.source_sha256=="652218003294aa3a582ace53c86f8c90749e643e9027d4b99a06cadb50d1d277" ? 29703 : size(model.coupling,1),"accuracy"=>errors,"synchronized_stage_seconds"=>stage_report,"prepared_EFV"=>result,"complete_EFV"=>complete,"independent_oracle"=>oracle_timing,"trajectory"=>sequence,"memory"=>memory))
            println(name," N=",n," ",T," median=",result["median_s"]," s ",errors); flush(stdout)
        end
    end
    report=Dict("schema"=>"PotentialFrameWork.GPUACEBenchmark","version"=>1,"options"=>opts,"julia"=>string(VERSION),
        "cpu"=>Sys.cpu_info()[1].model,"cpu_target"=>Sys.CPU_NAME,"gpu"=>storage===Array ? "omitted: CPU backend selected" : string(AMDGPU.device()),"hip_runtime"=>storage===Array ? "unavailable" : string(AMDGPU.HIP.runtime_version()),"os"=>string(Sys.KERNEL),"threads"=>Threads.nthreads(),"blas_threads"=>BLAS.get_num_threads(),
        "commit"=>strip(read(`git rev-parse HEAD`,String)),"dirty"=>!isempty(read(`git status --porcelain`,String)),"cases"=>records,"omitted"=>omitted,
        "versions"=>Dict(id.name=>string(pkgversion(mod)) for (id,mod) in Base.loaded_modules if id.name in ("AMDGPU","KernelAbstractions","AcceleratedKernels","SIMD","SpheriCart","BenchmarkTools")),
        "source_files_sha256"=>source_provenance(normpath(joinpath(@__DIR__,"../.."))),"kernel_release"=>Sys.isunix() ? strip(read(`uname -r`,String)) : "unavailable",
        "elapsed_seconds"=>time()-start,"boundary"=>"Prepared numeric EFV; synchronized; resident model/input; fresh cells at skin=0; owning=false GPU outputs borrow executor; CPU precision casting excluded", 
        "limitations"=>["Public legacy export verification required separately", "Oracle timing is in separate serialized fixtures; thread/SIMD sweeps use audit.jl and separate processes", "Device allocations not measured; host allocations are separate from buffer payloads"])
    open(opts["output"],"w") do io; TOML.print(io,report;sorted=true); end
end
if options(ARGS)["backend"]=="amd"
    @eval using AMDGPU
end
main()
