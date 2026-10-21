# No BenchmarkTools dependency required. Each measured call resets all outputs.
using PotentialFrameWork, Unitful, Statistics, InteractiveUtils, SIMD
include(joinpath(@__DIR__,"..","test","fixtures.jl"))
const NL=PotentialFrameWork.NeighbourList
const PP=PotentialFrameWork.PairPotentials
const REPEATS=parse(Int,get(ENV,"PFW_REPEATS","7"))

function geometry_scan(w::Workspace{W}) where W
    blocks=0; active=0; structural=0; candidates=0
    work=zeros(Int,length(w.keys))
    for task in eachindex(w.ranges), i in w.ranges[task], cell in w.adjacency[task]
        r=w.ranges[cell]
        for j=max(first(r),i+1):W:last(r)
            _,b=NL.geometry(w,i,j,last(r))
            blocks+=1; candidates+=W
            structural+=sum(b.structural); active+=sum(b.active)
            work[task]+=sum(b.active)
        end
    end
    (;blocks,active,structural,candidates,work)
end
function describe_samples(samples)
    (minimum(samples),median(samples),maximum(samples))
end
function run_case(n,cluster,mixed,width,parallel)
    s=fixture(n;L=max(4.1,cbrt(n/0.25)),cluster,mixed)
    p=potential(); w=Workspace(;width)
    for _=1:3; evaluate!(w,s,p;parallel); end
    stages=Dict(k=>Float64[] for k in keys(w.timings))
    samples=Float64[]; allocations=Int[]
    for _=1:REPEATS
        t=@timed evaluate!(w,s,p;parallel)
        push!(samples,t.time); push!(allocations,t.bytes)
        for (k,v) in w.timings; push!(stages[k],v); end
    end
    bound=bind_species(p,w.species_map)
    # Fixed-topology numeric reevaluation only. Public timings above always rebuild.
    for _=1:2; PP.execute!(w,bound,EnergyForcesVirial();parallel); end
    kernel=[(@timed PP.execute!(w,bound,EnergyForcesVirial();parallel)).time for _=1:REPEATS]
    kernelbytes=@allocated PP.execute!(w,bound,EnergyForcesVirial();parallel)
    geometry_scan(w)
    geometry=[(@timed geometry_scan(w)).time for _=1:REPEATS]
    occ=geometry_scan(w)
    total=describe_samples(samples).*1e3
    println("case N=$n cluster=$cluster mixed=$mixed W=$width parallel=$parallel")
    println("  public_ms(min,median,max)=",total," allocated_bytes_median=",median(allocations))
    println("  stages_ms_median=",sort([k=>median(v)*1e3 for (k,v) in stages];by=first))
    println("  fixed_topology_ms=",describe_samples(kernel).*1e3," allocated_bytes=",kernelbytes)
    println("  diagnostic_geometry_only_ms=",describe_samples(geometry).*1e3)
    println("  force_bytes=",sizeof(w.force)," workspace_bytes=",Base.summarysize(w),
        " task_partial_bytes=",sizeof(w.energies)+sizeof(w.virials))
    println("  phases=",length(w.phases)," tasks_per_phase=",length.(w.phases))
    println("  work_per_phase=",[sum(occ.work[phase]) for phase in w.phases],
        " task_work(min,max)=",extrema(occ.work))
    println("  blocks=",occ.blocks," active_pairs=",occ.active," structural_pairs=",occ.structural,
        " active_lane_fraction=",occ.active/max(1,occ.candidates),
        " structural_lane_fraction=",occ.structural/max(1,occ.candidates))
    @assert validate_schedule(w)
end

println("Julia=",VERSION," threads=",Threads.nthreads()," CPU=",Sys.cpu_info()[1].model," compiler_target=",Sys.CPU_NAME,
    " precision=Float64 repeats=",REPEATS," affinity=uncontrolled")
println("dependencies: AtomsBase=",pkgversion(AtomsBase)," SIMD=",pkgversion(SIMD)," Unitful=",pkgversion(Unitful))
println("Stage evaluation includes current-geometry packing, callback and immediate accumulation.")
println("Geometry-only scan is a separate diagnostic (also counts occupancy); not additive to stage timings.")
for n in (128,2048), (cluster,mixed) in ((false,false),(false,true),(true,true)), width in (1,4,8), parallel in (false,true)
    run_case(n,cluster,mixed,width,parallel)
end
if get(ENV,"PFW_LARGE","0")=="1"
    for width in (1,4,8), parallel in (false,true)
        run_case(8192,false,true,width,parallel)
    end
end
