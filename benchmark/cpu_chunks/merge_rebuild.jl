# Run from the repository root. ARGS: merge|rebuild OUTPUT.toml [samples=10] [grid=12]
include("threading.jl")

function observations(f, g, frames)
    functions = (f, g)
    records = (Dict{String,Any}[], Dict{String,Any}[])
    for run in -1:SAMPLES
        frame = frames[mod1(run, length(frames))]
        for index in (isodd(run) ? (1,2) : (2,1))
            call = functions[index]
            measured = @timed call(frame)
            run > 0 && push!(records[index], Dict("frame"=>mod1(run,length(frames)),
                "seconds"=>measured.time, "bytes"=>measured.bytes, "gc_s"=>measured.gctime))
        end
    end
    records
end

function active_graph(n)
    env = n.environment
    [sort!([(env.neighbours[k],Tuple(env.shifts[k])) for k in env.offsets[1,i]:env.offsets[end,i]-1
             if env.active[k]]) for i in eachindex(env.species)]
end

function rebuild_probes(sys, calc, ps, state, opts)

    stats = NL.ChunkTimings()
    public_call = @timed AC.calculate(EFV,sys,calc,ps,state; opts...,timings=stats)
    n = state.chunks.neighbours
    if hasproperty(stats, :rebuild)
        return Dict("environment_s"=>copy(n.environment.timings), "rebuild"=>copy(stats.rebuild),
                    "instrumented_public_s"=>public_call.time)
    end
    # Starting-revision component probes are independent of the timed public call.
    # Environment keys include nested extraction/sort/topology and snapshot traversal.
    input = A.ace_inputs(sys,calc.model;environment=A.FullEnvironment(;width=4),precision=Float64)
    order = @timed NL.chunk_neighbours(sys,input.environment,input.positions,input.cell)
    layout = @timed NL.form_numa_workers(order.value,1)
    routes = @timed NL.chunk_thread_plan(layout.value,opts.workers)
    metric(t) = Dict("seconds"=>t.time,"bytes"=>t.bytes,"gc_s"=>t.gctime)
    Dict("environment_s"=>copy(input.environment.timings),"instrumented_public_s"=>public_call.time,
         "independent"=>Dict("cell_order"=>metric(order),"layout"=>metric(layout),"routes"=>metric(routes)))
end

function run_campaign_case!(report, name, model_path, d)
    model, ps, st = A.load_resolved(model_path)
    displacement = [sin(i*(j+.31)) for j=1:3,i=1:size(d.positions,2)]
    frames = [TestSystem(d.positions .+ a .* displacement,d.cell,(true,true,true),
                         ChemicalSpecies.(d.species)) for a in (0.,1e-5,-1e-5)]
    calc = A.ACECalculator(model,ps,st;width=4,batch_size=1,slots=NTHREADS)
    state = A.ACEEvaluationState(calc)
    oracle = [AC.calculate(EFV,frame,calc;executor=:colored,parallel=false) for frame in frames]
    opts = (;executor=:threaded, workers=NTHREADS, groups=1)
    fresh(frame; kwargs...) = AC.calculate(EFV,frame,calc,ps,state;opts...,kwargs...)
    refresh(frame; kwargs...) = fresh(frame;reuse_topology=true,kwargs...)
    first_result = fresh(first(frames))
    reference_graph = active_graph(state.chunks.neighbours)
    permutation = copy(state.chunks.neighbours.p)
    for frame in frames
        fresh(frame)
        @assert active_graph(state.chunks.neighbours) == reference_graph
        @assert state.chunks.neighbours.p == permutation
    end
    record = Dict{String,Any}("name"=>name,"atoms"=>length(first(frames)),"model"=>relpath(model_path,ROOT),
        "model_sha256"=>sha(model_path),"width"=>4,"batch_size"=>1,"groups"=>1,
        "scatter_capacity"=>4096,"precision"=>"Float64","active_edges"=>sum(length,reference_graph))
    try
        if MODE == "merge"
            workers(frame; kwargs...) = refresh(frame; kwargs...)
            inline(frame; kwargs...) = refresh(frame; kwargs...,_inline_merge=true)
            plan = state.chunks.plan
            worker_results = workers.(frames)
            inline_results = inline.(frames)
            @assert state.chunks.plan === plan
            for (worker, local_merge, expected) in zip(worker_results, inline_results, oracle)
                check_result(worker,expected)
                @assert (worker.energy,worker.forces,worker.virial) ==
                        (local_merge.energy,local_merge.forces,local_merge.virial)
            end
            @assert all(name -> all(r -> r[name] != first(worker_results)[name],worker_results[2:end]),
                        (:energy,:forces,:virial))
            record["workers"], record["inline"] = observations(workers,inline,frames)
            for (policy,call) in (("workers",workers),("inline",inline))
                record[policy*"_diagnostics"] = [begin
                    stats = NL.ChunkTimings()
                    measured = @timed call(frame;timings=stats)
                    merge(timing_summary(stats,state.chunks,measured),
                          Dict("merge_waves"=>stats.merge_waves))
                end for frame in frames for _ in 1:4][2:end]
            end
        else
            for (frame, expected) in zip(frames,oracle)
                check_result(fresh(frame),expected)
                check_result(refresh(frame),expected)
            end
            record["fresh"], record["refresh"] = observations(fresh,refresh,frames)
            rebuild_probes(first(frames),calc,ps,state,opts)
            record["probes"] = [rebuild_probes(frame,calc,ps,state,opts) for frame in frames for _ in 1:3]
            moved = copy(d.positions)
            moved[:,1] .+= (1.01,0.79,-0.53)
            changed = TestSystem(moved,d.cell,(true,true,true),ChemicalSpecies.(d.species))
            for executor in (:chunks,:threaded), frame in (changed,first(frames))
                actual = AC.calculate(EFV,frame,calc,ps,state;executor,workers=NTHREADS)
                expected = AC.calculate(EFV,frame,calc;executor=:colored,parallel=false)
                check_result(actual,expected)
                @assert (active_graph(state.chunks.neighbours) == reference_graph) == (frame !== changed)
            end
            @assert first_result.state === state
            check_result(first_result,first(oracle))
            record["changed_graph_verified"] = true
        end
    finally
        close(state)
    end
    push!(report["cases"],record)
    save(report)
end

function campaign()
    @assert MODE in ("merge","rebuild") && SAMPLES>0 && GRID>=4
    BLAS.set_num_threads(1)
    affinity = strip(read(`taskset -pc $(getpid())`,String))
    topology = read(`lscpu -p=CPU,CORE,SOCKET,CACHE`,String)
    cpu_list = last(split(affinity,": "))
    cpus = reduce(vcat,map(split(cpu_list,",")) do part
        limits = parse.(Int,split(part,"-"))
        collect(first(limits):last(limits))
    end)
    rows = [split(row,",") for row in split(topology,'\n') if !isempty(row) && !startswith(row,"#")]
    physical = Set((row[2],row[3]) for row in rows if parse(Int,row[1]) in cpus)
    @assert length(cpus) == length(physical) == NTHREADS && Threads.nthreads(:interactive)==0
    report = Dict("mode"=>MODE,"cpu"=>Sys.cpu_info()[1].model,"julia"=>string(VERSION),
        "threads"=>NTHREADS,"blas_threads"=>BLAS.get_num_threads(),"affinity"=>affinity,
        "cpu_topology"=>topology,"governors"=>Dict(string(cpu)=>readchomp(path) for cpu in cpus
            for path in ("/sys/devices/system/cpu/cpu$cpu/cpufreq/scaling_governor",) if isfile(path)),
        "samples"=>SAMPLES,"grid"=>GRID,"source_sha256"=>source_hashes(),
        "starting_revision"=>"3ddd74c8fb90177d31414915ddff285e81c228ae",
        "harness_sha256"=>sha(@__FILE__),"threading_harness_sha256"=>sha(joinpath(@__DIR__,"threading.jl")),
        "fixture_sha256"=>sha(joinpath(ROOT,"test/fixtures.jl")),
        "workloads_sha256"=>sha(joinpath(ROOT,"benchmark/gpu_ace/workloads.jl")),
        "versions"=>Dict(id.name=>string(pkgversion(mod)) for (id,mod) in Base.loaded_modules
            if id.name in ("AtomsBase","AtomsCalculators","SIMD","SpheriCart","StaticArrays","Unitful",
                           "AcceleratedKernels","KernelAbstractions")),
        "cases"=>Dict{String,Any}[])
    models = joinpath(ROOT,"benchmark/cpu_ace/results/2026-11-10/models")
    run_campaign_case!(report,"Ti bcc",joinpath(models,"ti-resolved.toml"),bcc_geometry(GRID,[22]))
    run_campaign_case!(report,"Ti–Al–V bcc",joinpath(models,"tialv-resolved.toml"),bcc_geometry(GRID,[13,22,23]))
end
campaign()
