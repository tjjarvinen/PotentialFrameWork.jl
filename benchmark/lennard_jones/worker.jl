# One fresh Julia process per thread count. No package installation in measurements.
using PotentialFrameWork, AtomsBase, AtomsSystems, EmpiricalPotentials, StaticArrays
using Statistics, Unitful, Folds, NeighbourLists, SIMD
import AtomsCalculators as AC
using Pkg
include("common.jl")
using .ArgonComparison
const EX=PotentialFrameWork.Execution
const NL=PotentialFrameWork.NeighbourList
const SAMPLES=parse(Int,get(ENV,"PFW_LJ_REPEATS","7"))
const WARMUPS=parse(Int,get(ENV,"PFW_LJ_WARMUP","3"))
const CELLS=parse.(Int,split(get(ENV,"PFW_LJ_CELLS","4,6,10,16"),','))
SAMPLES>=3 && WARMUPS>=1 || error("use at least three samples and one warmup")
all(>=(4),CELLS) || error("periodic FCC systems require at least four cells/axis")

# Fixed-topology diagnostic, with the same boundary representation as AC calls.
# Internal Execution entrypoints are intentional and are NOT Verlet reuse APIs.
function fixed_call(w,bound,parallel)
    t=time_ns()
    EX.execute!(w,bound,EnergyForcesVirial();parallel)
    w.timings[:evaluation]=(time_ns()-t)*1e-9
    r=EX.outputs(w,EnergyForcesVirial())
    f=[SVector(r.forces[1,i],r.forces[2,i],r.forces[3,i]) for i in axes(r.forces,2)]
    (;energy=r.energy,forces=f,virial=SMatrix{3,3}(r.virial))
end

const STAGES=(:extraction,:sorting_packing,:topology,:schedule,:parameters,
              :evaluation,:partial_merge,:output_permutation)
function measure(io,method,regime,m,width,parallel,f,reference;workspace=nothing,setup_ms=0.)
    err=assert_agreement(f(),reference)
    for _=1:WARMUPS; f(); end
    GC.gc()
    times=Float64[]; bytes=Int[]; gc=Float64[]
    stages=Dict(k=>Float64[] for k in STAGES)
    for _=1:SAMPLES
        t=@timed f()
        push!(times,1e3*t.time); push!(bytes,t.bytes); push!(gc,1e3*t.gctime)
        if !isnothing(workspace)
            for key in STAGES
                used=regime=="rebuild" || key in (:evaluation,:partial_merge,:output_permutation)
                push!(stages[key],used ? workspace.timings[key]*1e3 : 0.)
            end
        end
    end
    assert_agreement(f(),reference) # reset/cached paths must remain numerically correct
    n=4*m^3
    force_bytes=isnothing(workspace) ? 0 : sizeof(workspace.force)
    phases=isnothing(workspace) ? 0 : length(workspace.phases)
    stage_medians=[isnothing(workspace) ? NaN : median(stages[k]) for k in STAGES]
    fields=[Threads.nthreads(),n,m,method,regime,width,parallel,SAMPLES,
            minimum(times),median(times),maximum(times),median(bytes),median(gc),
            err.energy,err.force,err.virial,setup_ms,force_bytes,phases,stage_medians...]
    println(io,join(fields,',')); flush(io)
    println("N=$n $method $regime W=$width parallel=$parallel median_ms=$(median(times)) allocated_bytes=$(median(bytes)) errors=$err")
    if !isnothing(workspace)
        @assert validate_schedule(workspace)
        println("  phases=$phases tasks_per_phase=$(length.(workspace.phases)) force_bytes=$force_bytes workspace_bytes=$(Base.summarysize(workspace))")
    end
    flush(stdout)
end

# Match the reference assembly's actual constructor on each supported backend.
function reference_list(sys,radius)
    if pkgversion(NeighbourLists)>=v"0.6"
        return NeighbourLists.PairList(sys,radius*u"Å";int_type=Int)
    end
    NeighbourLists.PairList(sys,radius*u"Å")
end

function main()
    output=isempty(ARGS) ? joinpath(@__DIR__,"results") : abspath(ARGS[1])
    mkpath(output)
    threads=Threads.nthreads()
    logfile=joinpath(output,"metadata-threads-$threads.txt")
    metadata="Julia=$(VERSION) default_threads=$threads CPU=$(Sys.cpu_info()[1].model) compiler_target=$(Sys.CPU_NAME) affinity=uncontrolled precision=Float64 warmups=$WARMUPS samples=$SAMPLES cells=$CELLS"
    println(metadata)
    open(logfile,"w") do io
        println(io,metadata)
        println(io,"epsilon=0.010323 eV sigma=3.405 Å rmin=2^(1/6)*sigma cutoff=2.5*sigma FCC_a=5.26 Å deterministic_rattle=0.04 Å")
        println(io,"All timed calls return energy, original-order SVector forces and SMatrix virial.")
        println(io,"rebuild includes topology each call; fixed excludes one-time topology/parameter setup (setup_ms), positions remain unchanged.")
        println(io,"EP SequentialEx ntasks=1; EP ThreadedEx ntasks=default_threads. NeighbourLists uses its default CPU build backend.")
        println(io,"EP site assembly uses full-size per-task force partials; PFW uses one shared sorted force workspace.")
        for entry in sort(collect(values(Pkg.dependencies()));by=x->x.name)
            println(io,entry.name,"=",entry.version)
        end
    end
    csv=joinpath(output,"threads-$threads.csv")
    open(csv,"w") do io
        columns=["threads","N","fcc_cells","implementation","regime","width","parallel","samples",
                 "min_ms","median_ms","max_ms","allocated_bytes_median","gc_ms_median",
                 "energy_abs_error_eV","force_max_abs_error_eV_per_A","virial_max_abs_error_eV",
                 "setup_ms_excluded","sorted_force_bytes","phase_count",[string(k)*"_median_ms" for k in STAGES]...]
        println(io,join(columns,','))
        ours,ep=argon_models()
        for m in CELLS
            sys=argon_fcc(m)
            ref=AC.energy_forces_virial(sys,ep;executor=SequentialEx(),ntasks=1)
            # Complete public calls include list creation in both packages.
            for parallel in (false,true)
                ex=parallel ? ThreadedEx() : SequentialEx()
                tasks=parallel ? threads : 1
                f=()->AC.energy_forces_virial(sys,ep;executor=ex,ntasks=tasks)
                measure(io,"EmpiricalPotentials","rebuild",m,0,parallel,f,ref)
            end
            for (width,parallel) in ((1,false),(4,false),(4,true),(8,false),(8,true))
                w=Workspace(;width)
                f=()->AC.energy_forces_virial(sys,ours;neighbourlist=w,parallel)
                measure(io,"PotentialFrameWork","rebuild",m,width,parallel,f,ref;workspace=w)
            end
            # Separately labelled unchanged-geometry diagnostics; no setup is hidden
            # in the complete public rows above.
            built=@timed reference_list(sys,ours.radius)
            nlist=built.value
            for parallel in (false,true)
                ex=parallel ? ThreadedEx() : SequentialEx()
                tasks=parallel ? threads : 1
                f=()->AC.energy_forces_virial(sys,ep;nlist,executor=ex,ntasks=tasks)
                measure(io,"EmpiricalPotentials","fixed",m,0,parallel,f,ref;setup_ms=1e3*built.time)
            end
            for (width,parallel) in ((1,false),(4,false),(4,true),(8,false),(8,true))
                w=Workspace(;width)
                # Compile/warm preparation first, then measure this one setup.
                NL.prepare!(w,sys,ours.radius)
                setup=@timed NL.prepare!(w,sys,ours.radius)
                bound=bind_species(ours,w.species_map)
                f=()->fixed_call(w,bound,parallel)
                measure(io,"PotentialFrameWork","fixed",m,width,parallel,f,ref;
                        workspace=w,setup_ms=1e3*setup.time)
            end
        end
    end
end
main()
