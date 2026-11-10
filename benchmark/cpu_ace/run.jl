# CPU-only BenchmarkTools worker; no GPU package/device import.
include("common.jl")
using Profile
function options(args)
    d=Dict("preset"=>"quick","data"=>joinpath(@__DIR__,"results/2026-11-10/models"),
        "models"=>"small-tial,medium-ti,large-tial","oracle-data"=>"","sizes"=>"","precisions"=>"Float64,Float32",
        "samples"=>"10","seconds"=>"3","output"=>"cpu-ace-results.toml","widths"=>"4",
        "diagnostics"=>"true","profile"=>"","codegen"=>"","parallel"=>"false","historical"=>"false","pace"=>get(ENV,"PFW_LAMMPS","/opt/lammps-static/bin/lmp"))
    iseven(length(args)) || error("use --option value pairs")
    for i=1:2:length(args)
        k=replace(args[i],r"^--"=>""); haskey(d,k) || error("unknown option $k"); d[k]=args[i+1]
    end
    d["preset"] in ("quick","full") || error("preset")
    d
end
function main()
    opts=options(ARGS); BLAS.set_num_threads(1)
    samples=parse(Int,opts["samples"]); seconds=parse(Float64,opts["seconds"])
    samples>0 && seconds>0 || error("positive samples/seconds required")
    sizes=parse.(Int,split(isempty(opts["sizes"]) ? (opts["preset"]=="quick" ? "4,18" : "4,8,18") : opts["sizes"],','))
    widths=parse.(Int,split(opts["widths"],',')); parallel=opts["parallel"]=="true"
    historical=opts["historical"]=="true"; records=Dict[]; omissions=Dict[]
    source=source_provenance(dirname(dirname(pathof(PotentialFrameWork))))
    # Persist partial results after every completed case, including validation failures.
    report=Dict("schema"=>"PFW.CPUACE.Stage1.v1","options"=>opts,"julia"=>string(VERSION),
        "cpu"=>Sys.cpu_info()[1].model,"llvm_version"=>string(Base.libllvm_version),"cpu_target"=>Sys.CPU_NAME,"threads"=>Threads.nthreads(),
        "blas_threads"=>BLAS.get_num_threads(),"affinity"=>"uncontrolled; inspect process environment",
        "source_files_sha256"=>source,
        "harness_sha256"=>Dict(f=>bytes2hex(sha256(read(joinpath(@__DIR__,f)))) for f in ("run.jl","common.jl","generate_modern.jl","oracle_modern.jl","legacy_oracle.jl","features.jl","suite.jl")),"loaded_package_source"=>pathof(PotentialFrameWork),
        "versions"=>Dict(id.name=>string(pkgversion(mod)) for (id,mod) in Base.loaded_modules if id.name in ("SIMD","SpheriCart","BenchmarkTools","StaticArrays","KernelAbstractions","AcceleratedKernels")),
        "cases"=>records,"omissions"=>omissions)
    output=abspath(opts["output"]); mkpath(dirname(output))
    save()=open(io->TOML.print(io,report;sorted=true),output,"w")
    for name in split(opts["models"],',')
        path=joinpath(opts["data"],name*"-resolved.toml"); oraclepath=joinpath(isempty(opts["oracle-data"]) ? opts["data"] : opts["oracle-data"],name*"-oracle.toml")
        if !isfile(path)
            push!(omissions,Dict("model"=>name,"status"=>"skipped","reason"=>"resolved model missing: $path")); continue
        end
        model_load_s=@elapsed model,ps,_=A.load_resolved(path)
        oracle=isfile(oraclepath) ? TOML.parsefile(oraclepath) : nothing
        for g in sizes, Tname in split(opts["precisions"],','), width in widths
            opts["preset"]=="quick" && name!="small-tial" && g!=first(sizes) && continue
            T=Tname=="Float64" ? Float64 : Tname=="Float32" ? Float32 : error("precision")
            adaptation_s=@elapsed b=A.adapt_model(model,ps;precision=T)
            d=bcc_geometry(g,model.atomic_numbers)
            sys=structure(d;T); env=A.FullEnvironment(;width)
            elapsed=@elapsed input=A.ace_inputs(sys,b.model;environment=env,precision=T)
            X=T.(input.positions); H=T.(input.cell)
            w=A.CPUExecutor(b.model,env;precision=T)
            algorithms=historical ? (:historical,) : (:reference,:optimized)
            ref=nothing
            for algorithm in algorithms
                kw=historical ? (;parallel) : (;parallel,product_algorithm=algorithm)
                f=()->A.execute_prepared!(w,X,H,b.parameters;kw...)
                first_prepared_s=@elapsed f()
                r=A.prepared_result(w)
                if isnothing(ref); ref=r; else; agree(r,ref,T); end
                record=Dict("model"=>name,"g"=>g,"atoms"=>length(sys),"precision"=>Tname,"width"=>width,
                    "algorithm"=>string(algorithm),"backend"=>"CPU",
                    "model_load_s"=>model_load_s,"adaptation_s"=>adaptation_s,"first_prepared_call_s"=>first_prepared_s,"metadata"=>metadata(b.model,b.parameters),
                    "model_file_sha256"=>bytes2hex(sha256(read(path))),
                    "geometry_sha256"=>bytes2hex(sha256(reinterpret(UInt8,vec(X)))),
                    "species_sha256"=>bytes2hex(sha256(string(d.species))),"cell_sha256"=>bytes2hex(sha256(reinterpret(UInt8,vec(d.cell)))),
                    "edges"=>count(env.active),"packed_lanes"=>length(env.active),"phase_sizes"=>length.(env.phases),
                    "scratch_payload_bytes"=>sum(s->sum(sizeof(getfield(s,k)) for k in fieldnames(typeof(s))),w.scratch),
                    "single_shared_force_bytes"=>sizeof(w.force),"prepared"=>measure(f,samples,seconds),
                    "reference_errors"=>accuracy(r,ref),"first_preparation_s"=>elapsed,"preparation_stage_s"=>Dict(string(k)=>v for (k,v) in env.timings))
                if !isnothing(oracle)
                    o=findfirst(x->x["g"]==g,oracle["records"])
                    if !isnothing(o)
                        raw=oracle["records"][o]
                        external=(;energy=raw["energy"],forces=reshape(raw["forces"],3,:),virial=reshape(raw["virial"],3,3))
                        agree(r,external,T)
                        record["external_Float64_errors"]=accuracy(r,external)
                        record["external_boundary"]=raw["boundary"]
                        record["external_evaluator"]="ACEpotentials "*oracle["ACEpotentials"]*(oracle["ACEpotentials"]=="0.6.6" ? " / ACE1 analytic" : " / modern")
                        @assert oracle["source_sha256"]==model.source_sha256
                        record["external_raw_s"]=raw["raw_ns"].*1e-9
                        record["external_precision_comparison"]=T===Float64 ? "matched Float64" : "application sensitivity against Float64; external Float32 unsupported"
                    end
                end
                calc=historical ? A.ACECalculator(b.model,b.parameters;precision=T,width) : A.ACECalculator(b.model,b.parameters;precision=T,width,product_algorithm=algorithm)
                complete=()->AC.energy_forces_virial(sys,calc;parallel)
                record["complete_public"]=measure(complete,samples,seconds)
                # Modestly perturbed replay, always fresh connectivity per call.
                frames=[structure((;positions=d.positions .+ a.*[sin(i*(j+.31)) for j=1:3,i=1:length(sys)],cell=d.cell,species=d.species);T) for a in (0.,.015,.03)]
                replay=()->foreach(frame->AC.energy_forces_virial(frame,calc;parallel),frames)
                record["repeated_frames"]=merge(measure(replay,samples,seconds),Dict("frames"=>3,"connectivity"=>"fresh rebuild every frame; includes public owning outputs"))
                # Same operation can be profiled separately, without stage timers.
                if opts["diagnostics"]=="true" && g==first(sizes) && !historical
                    A.prepare_executor!(w,env)
                    codepath=isempty(opts["codegen"]) ? "" : joinpath(opts["codegen"],name*"-"*Tname*"-"*string(algorithm)*"-w"*string(width))
                    record["isolated_site_stages"]=diagnostics(w,X,H,b.parameters,algorithm,samples,seconds,codepath)
                    prepare=()->begin
                        inp=A.ace_inputs(sys,b.model;environment=env,precision=T)
                        A.prepare_executor!(w,env)
                    end
                    record["preparation"]=measure(prepare,samples,seconds)
                    record["executor_preparation"]=measure(()->A.prepare_executor!(w,env),samples,seconds)
                    record["preparation_stage_s"]=Dict(string(k)=>v for (k,v) in env.timings)
                    stage_observations=Dict[]
                    for _=1:samples
                        observation=@timed prepare()
                        push!(stage_observations,Dict("seconds"=>observation.time,"bytes"=>observation.bytes,
                            "stages_s"=>Dict(string(k)=>v for (k,v) in env.timings)))
                    end
                    record["preparation_diagnostics_raw"]=stage_observations
                    extract=()->begin
                        positions=[U.ustrip(U.u"Å",AB.position(sys,i)[a]) for a=1:3,i=1:length(sys)]
                        cell=[U.ustrip(U.u"Å",AB.cell_vectors(sys)[j][a]) for a=1:3,j=1:3]
                        (;positions,cell)
                    end
                    record["boundary_extraction"]=measure(extract,samples,seconds)
                    record["neighbour_construction"]=measure(()->NL.neighbour_list(sys,A.NI.search_cutoff(b.model);width,precision=T),samples,seconds)
                    if !isempty(opts["profile"]) && algorithm===:optimized
                        Profile.clear()
                        Profile.@profile for _=1:20; f(); end
                        mkpath(opts["profile"])
                        open(joinpath(opts["profile"],name*"-"*Tname*"-w"*string(width)*".txt"),"w") do io
                            Profile.print(io;format=:flat,sortedby=:count)
                        end
                    end
                end
                push!(records,record); save()
                println(name," ",length(sys)," ",Tname," ",algorithm," prepared ms=",record["prepared"]["median_s"]*1000); flush(stdout)
            end
        end
    end
    # Optional PACE probe only. No timing is claimed before spline equivalence.
    pace=opts["pace"]
    reason=if isfile(pace)
        text=try read(`$pace -h`,String) catch err; sprint(showerror,err) end
        report["pace_probe"]=text
        report["pace_executable_sha256"]=bytes2hex(sha256(read(pace)))
        "Executable probed; equivalent exported spline interpolation and pair tables not validated in stage 1. No PACE comparison launched."
    else
        "LAMMPS executable unavailable at $pace"
    end
    push!(omissions,Dict("model"=>"PACE","status"=>"unsupported","reason"=>reason))
    save()
end
main()
