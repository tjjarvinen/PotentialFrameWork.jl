# Run from the repository root with Julia on one logical CPU per physical core.
# ARGS: baseline|updated|diagnostics OUTPUT.toml [samples=7] [grid=12]
using PotentialFrameWork, AtomsBase, SpheriCart, Unitful, StaticArrays, LinearAlgebra, Statistics, TOML, SHA
import AtomsCalculators as AC
include("../../test/fixtures.jl")
include("../gpu_ace/workloads.jl")
const A = PotentialFrameWork.ACE
const NL = PotentialFrameWork.NeighbourList
const ROOT = normpath(joinpath(@__DIR__, "../.."))
const MODE, OUTPUT = ARGS[1:2]
const SAMPLES = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 7
const GRID = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 12
const NTHREADS = Threads.nthreads(:default)
const EFV = (AC.Energy(), AC.Forces(), AC.Virial())

function sample(f)
    f(); f()
    GC.gc()
    runs = [@timed f() for _ in 1:SAMPLES]
    times = getproperty.(runs, :time)
    Dict("seconds"=>times, "median_s"=>median(times), "min_s"=>minimum(times),
         "max_s"=>maximum(times), "bytes"=>getproperty.(runs, :bytes),
         "gc_s"=>getproperty.(runs, :gctime))
end

function replay(f, frames)
    i = Ref(0)
    () -> begin
        i[] = mod1(i[]+1, length(frames))
        f(frames[i[]])
    end
end

function sha(path)
    open(io -> bytes2hex(sha256(io)), path)
end

function source_hashes()
    root = dirname(dirname(pathof(PotentialFrameWork)))
    Dict(relpath(path, root) => sha(path) for dir in ("src","ext")
         for (base, _, files) in walkdir(joinpath(root,dir)) for file in files
         for path in (joinpath(base,file),) if endswith(file,".jl"))
end

function timing_summary(stats, w, measured)
    n = w.neighbours
    env = n.environment
    edges = [sum(group) do k
        i = n.p[k]
        count(view(env.active, env.offsets[1,i]:env.offsets[end,i]-1))
    end for group in w.groups]
    plan = length(stats.site)>1 ? w.plan : nothing
    slots = isnothing(plan) ? (isempty(w.chunks) ? 0 : 1) : plan.workers * plan.batches
    Dict("public_wall_s"=>measured.time, "public_allocated_bytes"=>measured.bytes,
         "pass_wall_s"=>stats.pass * 1e-9, "site_barriers_s"=>stats.site_barrier * 1e-9,
         "merge_barriers_s"=>stats.merge_barrier * 1e-9, "export_s"=>stats.export_time * 1e-9,
         "site_slot_s"=>stats.site .* 1e-9, "merge_slot_s"=>stats.merge .* 1e-9,
         "site_slot_noncompute_s"=>[max(0, Float64(stats.site_barrier)-t) * 1e-9 for t in stats.site],
         "slot_ack_wait_s"=>stats.slot_wait .* 1e-9, "allocated_slots"=>slots,
         "peak_active_slots"=>stats.peak_slots, "slot_capacity"=>stats.capacity,
         "force_scatter_bytes"=>(haskey(w.results,:forces) ? slots * stats.capacity * sizeof(eltype(w.results.forces)) : 0),
         "wave_count"=>isnothing(plan) ? 0 : length(plan.waves),
         "chunks"=>length(w.chunks), "group_sites"=>length.(w.groups), "group_active_edges"=>edges)
end

function colored_prepared(sys, calc, props)
    X,H = A.ace_geometry(sys)
    A.execute_prepared!(calc.executor,X,H,calc.parameters;properties=props,parallel=true)
    r = A.prepared_result(calc.executor)
    (;energy=r.energy*u"eV",forces=A.force_vectors(r.forces*u"eV/Å"),virial=A.static_virial(r.virial*u"eV"))
end

function check_result(x, y)
    for (name, unit) in ((:energy,u"eV"),(:forces,u"eV/Å"),(:virial,u"eV"))
        haskey(y,name) || continue
        numeric(v::Number) = ustrip(unit,v)
        numeric(v::AbstractArray) = map(numeric,v)
        @assert isapprox(numeric(x[name]),numeric(y[name]); atol=2e-8, rtol=2e-8)
    end
end

function run_case!(report, name, model_path, d, periodic, small)
    small && NTHREADS > 4 && return
    model, ps, st = A.load_resolved(model_path)
    frames = [TestSystem(d.positions .+ a .* [sin(i*(j+.31)) for j=1:3,i=1:size(d.positions,2)],
                         d.cell,periodic,ChemicalSpecies.(d.species)) for a in (0.,1e-5,-1e-5)]
    calc = A.ACECalculator(model,ps,st; width=4, batch_size=1, slots=NTHREADS)
    serial = A.ACEEvaluationState(calc)
    input = A.ace_inputs(first(frames), calc.model; environment=A.FullEnvironment(;width=4), precision=Float64)
    env = input.environment
    # Check replay geometry keeps the exact graph, image shifts and atom mapping.
    baseline = NL.chunk_neighbours(first(frames),env,input.positions,input.cell)
    for sys in frames[2:end]
        i = A.ace_inputs(sys,calc.model;environment=A.FullEnvironment(;width=4),precision=Float64)
        n = NL.chunk_neighbours(sys,i.environment,i.positions,i.cell)
        @assert n.p == baseline.p && n.ranges == baseline.ranges
        @assert i.environment.offsets == env.offsets && i.environment.neighbours == env.neighbours
        @assert i.environment.shifts == env.shifts && i.environment.active == env.active
    end
    record = Dict{String,Any}("name"=>name,"atoms"=>length(first(frames)),"active_edges"=>count(env.active),
        "packed_lanes"=>length(env.active),"occupied_cells"=>length(baseline.ranges),"phase_sizes"=>length.(env.phases),
        "model"=>relpath(model_path,ROOT),"model_sha256"=>sha(model_path),"width"=>4,"batch_size"=>1,
        "scatter_capacity"=>4096,"precision"=>"Float64","runs"=>Dict{String,Any}[])
    push!(report["cases"],record)
    save(report)
    if MODE != "diagnostics"
        record["cell_order_rebuild"] = sample(() -> NL.chunk_neighbours(first(frames),env,input.positions,input.cell))
        record["environment_rebuild"] = sample(replay(frames) do sys
            A.ace_inputs(sys,calc.model;environment=serial.environment,precision=Float64)
        end)
        record["chunk_layout"] = sample(() -> NL.form_numa_workers(baseline,1))
    end

    jobs = small ? ((:energy, (AC.Energy(),)), (:forces, (AC.Forces(),))) : ((:efv,EFV),)
    for (request, props) in jobs
        oracle = [AC.calculate(props,s,calc;executor=:colored,parallel=false) for s in frames]
        function run(method, policy, groups, f)
            replayed = f.(frames)
            for (x, y) in zip(replayed,oracle); check_result(x,y); end
            for name in keys(first(oracle))
                name === :state && continue
                @assert any(x -> x[name] != first(replayed)[name], replayed[2:end])
            end
            timings = sample(replay(f,frames))
            push!(record["runs"],merge(timings,Dict("method"=>method,"policy"=>policy,
                "groups"=>groups,"request"=>string(request))))
            save(report)
        end
        if MODE == "baseline"
            run("serial","rebuild",1,s -> AC.calculate(props,s,calc,ps,serial))
            run("colored","rebuild",1,s -> AC.calculate(props,s,calc;executor=:colored,parallel=true))
        elseif MODE == "updated"
            calc.executor.environment=env
            A.prepare_executor!(calc.executor,env)
            run("colored","prepared reuse",1,s -> colored_prepared(s,calc,props))
        end
        if MODE != "baseline"
            for groups in (1, 4)
                state = A.ACEEvaluationState(calc)
                options = (;executor=:threaded,workers=NTHREADS,groups)
                try
                    if MODE == "updated"
                        run("threaded","rebuild",groups,s -> AC.calculate(props,s,calc,ps,state;options...))
                        AC.calculate(props,first(frames),calc,ps,state;options...)
                        run("threaded","reuse",groups,s -> AC.calculate(props,s,calc,ps,state;options...,reuse_topology=true))
                    else
                        AC.calculate(props,first(frames),calc,ps,state;options...)
                        for reuse in (true,false)
                            AC.calculate(props,frames[2],calc,ps,state;options...,reuse_topology=reuse,timings=NL.ChunkTimings())
                            stats = NL.ChunkTimings()
                            measured = @timed AC.calculate(props,frames[3],calc,ps,state;options...,
                                                          reuse_topology=reuse,timings=stats)
                            push!(record["runs"],merge(timing_summary(stats,state.chunks,measured),Dict(
                                "method"=>"threaded","groups"=>groups,"request"=>string(request),
                                "policy"=>reuse ? "reuse" : "rebuild")))
                            save(report)
                        end
                    end
                finally
                    close(state)
                end
            end
        end
    end
end

save(report) = open(io -> TOML.print(io,report;sorted=true),OUTPUT,"w")

function main()
    @assert MODE in ("baseline","updated","diagnostics") && SAMPLES>0 && GRID>=4
    BLAS.set_num_threads(1)
    report = Dict("mode"=>MODE,"cpu"=>Sys.cpu_info()[1].model,"julia"=>string(VERSION),
        "threads"=>NTHREADS,"interactive_threads"=>Threads.nthreads(:interactive),"blas_threads"=>BLAS.get_num_threads(),
        "affinity"=>strip(read(`taskset -pc $(getpid())`,String)),"samples"=>SAMPLES,"grid"=>GRID,
        "baseline_commit"=>"90d726d08c21aa3e7bbd9b60dc0d17d31369c3f5","source_sha256"=>source_hashes(),
        "llvm"=>string(Base.libllvm_version),"cpu_target"=>Sys.CPU_NAME,
        "fixture_sha256"=>sha(joinpath(@__DIR__,"../../test/fixtures.jl")),
        "harness_sha256"=>sha(@__FILE__), "workloads_sha256"=>sha(joinpath(@__DIR__,"../gpu_ace/workloads.jl")),
        "versions"=>Dict(id.name=>string(pkgversion(mod)) for (id,mod) in Base.loaded_modules
            if id.name in ("AtomsBase","AtomsCalculators","SIMD","SpheriCart","StaticArrays","Unitful","AcceleratedKernels","KernelAbstractions")),
        "cases"=>Dict{String,Any}[])
    models = joinpath(ROOT,"benchmark/cpu_ace/results/2026-11-10/models")
    ti = joinpath(models,"ti-resolved.toml")
    alloy = joinpath(models,"tialv-resolved.toml")
    dimer = (positions=[0. 2.6;0. 0.1;0. 0.],cell=20 .* Matrix{Float64}(I,3,3),species=[22,22])
    run_case!(report,"Ti2 dimer",ti,dimer,(false,false,false),true)
    run_case!(report,"Ti 16-atom periodic",ti,bcc_geometry(2,[22]),(true,true,true),true)
    run_case!(report,"Ti bcc",ti,bcc_geometry(GRID,[22]),(true,true,true),false)
    run_case!(report,"Ti–Al–V bcc",alloy,bcc_geometry(GRID,[13,22,23]),(true,true,true),false)
end
main()
