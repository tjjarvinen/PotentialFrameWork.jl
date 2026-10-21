"""User-facing comparison of neighbour-search strategies without a potential."""
module NeighbourBenchmarks
using AtomsBase: AbstractSystem
using Unitful
using ..NeighbourList: Workspace, prepare!, lane
using ..Execution: geometry
export benchmark_neighbourlist

@inline function task_pairs(w::Workspace{W}, task) where W
    pairs=0
    for i in w.ranges[task], image in w.images[task]
        r=w.ranges[image.cell]
        for j=max(first(r),i+(image.self ? 0 : 1)):W:last(r)
            _,block=geometry(w,i,j,last(r),image)
            for k=1:W
                pairs+=lane(block.active,k)
            end
        end
    end
    pairs
end

function scan!(counts,w,parallel)
    resize!(counts,length(w.ranges))
    if parallel
        for phase in w.phases
            Threads.@threads :dynamic for k in eachindex(phase)
                task=phase[k]
                counts[task]=task_pairs(w,task)
            end
        end
    else
        for task in eachindex(w.ranges)
            counts[task]=task_pairs(w,task)
        end
    end
    sum(counts)
end

function geometry_pass!(w,counts,sys,cutoff,parallel)
    start=time_ns()
    prepare!(w,sys,cutoff)
    prepared=time_ns()
    pairs=scan!(counts,w,parallel)
    finished=time_ns()
    (;pairs,preparation=(prepared-start)*1e-9,traversal=(finished-prepared)*1e-9)
end

measure!(w,counts,sys,cutoff,parallel)=@timed geometry_pass!(w,counts,sys,cutoff,parallel)

# Avoid adding a Statistics/BenchmarkTools dependency for a small diagnostic.
function median_value(values)
    sorted=sort(values)
    n=length(sorted)
    isodd(n) ? sorted[(n+1)÷2] : (sorted[n÷2]+sorted[n÷2+1])/2
end
function summarize(samples)
    seconds=[s.time for s in samples]
    bytes=[s.bytes for s in samples]
    (;minimum=minimum(seconds),median=median_value(seconds),maximum=maximum(seconds),
      preparation_median=median_value([s.value.preparation for s in samples]),
      traversal_median=median_value([s.value.traversal for s in samples]),
      bytes_median=median_value(bytes),seconds,bytes)
end

"""
    benchmark_neighbourlist(sys, cutoff; width=4, precision=Float64,
                           parallel=false, samples=20, warmup=3, io=stdout)

Compare forced brute-force and cell-list paths on an AtomsBase system. `cutoff`
is a numeric length in Å or a unitful length. Each sample rebuilds topology and
schedule, then traverses all candidate SIMD blocks and counts active pair images.
Separate private workspaces are reused; the caller's system/settings are unchanged.
Warmups exclude compilation and grow reusable buffers before measurement. Sample
order alternates to reduce order bias. GC remains enabled; times are in seconds,
and allocations count Julia-visible bytes, including preparation and traversal.

Print min/median/max total time, median preparation/traversal times and allocated
bytes. Return `(;bruteforce,cells,preferred,cells_over_bruteforce,pairs,metadata)`;
each method includes raw `seconds` and `bytes` samples. Set `io=nothing` to suppress
printing. `preferred` is the method with lower median geometry time (or
`:equivalent` for an empty system); the ratio is cell time / brute-force time.

This is a GEOMETRY benchmark: it excludes potential callbacks, force/virial
accumulation and output permutation. Benchmark the real potential separately
when those costs matter. Results apply to this system, cutoff, width, precision,
thread setting and hardware; they do not establish a universal atom threshold.
Brute force has one serial task even with `parallel=true`. Empty systems have
identical empty topology for both settings. Run without other benchmarks/tests.
"""
function benchmark_neighbourlist(sys::AbstractSystem{3},cutoff;
        width=4,precision=Float64,parallel=false,samples=20,warmup=3,io=stdout)
    samples isa Integer && samples>0 || throw(ArgumentError("samples must be a positive integer"))
    warmup isa Integer && warmup>0 || throw(ArgumentError("warmup must be a positive integer"))
    c=cutoff isa Unitful.Length ? ustrip(u"Å",cutoff) : cutoff
    c isa Real || throw(ArgumentError("cutoff must be numeric Å or a unitful length"))
    isfinite(c) && c>0 || throw(ArgumentError("cutoff must be finite and positive"))
    n=length(sys)
    brute=Workspace(;width,precision,bruteforce_threshold=n)
    cells=Workspace(;width,precision,bruteforce_threshold=0)
    workspaces=(brute,cells)
    counts=(Int[],Int[])
    for _=1:warmup, method=1:2
        measure!(workspaces[method],counts[method],sys,c,parallel)
    end
    expected=sum(counts[1])
    expected==sum(counts[2]) || error("neighbour methods disagree on active pair-image count")
    measurements=ntuple(_ -> Vector{NamedTuple}(undef,samples),2)
    for sample=1:samples
        order=isodd(sample) ? (1,2) : (2,1)
        for method in order
            result=measure!(workspaces[method],counts[method],sys,c,parallel)
            result.value.pairs==expected || error("active pair-image count changed during benchmark")
            measurements[method][sample]=result
        end
    end
    b=summarize(measurements[1]); normal=summarize(measurements[2])
    preferred=n==0 ? :equivalent : b.median<=normal.median ? :bruteforce : :cells
    metadata=(;atoms=n,cutoff_angstrom=precision(c),width,precision,parallel,
              threads=Threads.nthreads(),julia_version=VERSION,cpu=Sys.cpu_info()[1].model,samples,warmup)
    report=(;bruteforce=b,cells=normal,preferred,
             cells_over_bruteforce=normal.median/b.median,pairs=expected,metadata)
    isnothing(io) || print_report(io,report)
    report
end

function print_report(io,r)
    m=r.metadata
    println(io,"Neighbour geometry benchmark: ",m.atoms," atoms, cutoff=",m.cutoff_angstrom,
        " Å, ",r.pairs," active pair images")
    println(io,"Julia ",m.julia_version,", ",m.cpu,", threads=",m.threads,
        ", parallel=",m.parallel,", width=",m.width,", precision=",m.precision)
    println(io,m.samples," samples per method after ",m.warmup," warmups; times in ms")
    println(io,"method       min / median / max; preparation median; traversal median; allocated bytes median")
    for (name,s) in (("bruteforce",r.bruteforce),("cells",r.cells))
        println(io,rpad(name,13),join(round.([s.minimum,s.median,s.maximum].*1000;digits=4)," / "),
            "; ",round(s.preparation_median*1000;digits=4),"; ",round(s.traversal_median*1000;digits=4),
            "; ",s.bytes_median)
    end
    println(io,"Lower median geometry time: ",r.preferred,"; cells/bruteforce = ",
        round(r.cells_over_bruteforce;digits=3))
    println(io,"Excludes potential evaluation, force/virial accumulation and output conversion.")
    if m.atoms>0
        println(io,"Select with Workspace(bruteforce_threshold=",m.atoms,
            ") for brute force, or Workspace(bruteforce_threshold=0) for cells.")
    end
end
end
