# Worker tasks own their inbox and completion channel. They hold neither the
# workgroup nor the last pass's callback while idle.
mutable struct ChunkWorkerPool
    inboxes::Vector{Channel{Any}}
    completed::Vector{Channel{Any}}
    tasks::Vector{Task}
    storage::Any
end

function chunk_worker_loop(inbox,completed,slot)
    for job in inbox
        failure=try
            # A later pass can submit a callback defined after pool startup.
            Base.invokelatest(job,slot)
            nothing
        catch exception
            CapturedException(exception,catch_backtrace())
        end
        job=nothing
        put!(completed,failure)
    end
end

function ChunkWorkerPool(workers)
    inboxes=[Channel{Any}(1) for _ in 1:workers]
    completed=[Channel{Any}(1) for _ in 1:workers]
    tasks=Task[]
    try
        for slot in 1:workers
            push!(tasks,Threads.@spawn :default chunk_worker_loop(inboxes[slot],completed[slot],slot))
        end
    catch
        foreach(close,inboxes)
        foreach(wait,tasks)
        rethrow()
    end
    ChunkWorkerPool(inboxes,completed,tasks,nothing)
end

function Base.close(pool::ChunkWorkerPool)
    foreach(close,pool.inboxes)
    foreach(wait,pool.tasks)
    pool.storage=nothing
    nothing
end

function dispatch_chunk_workers!(job,pool::ChunkWorkerPool,workers)
    dispatched=0
    failure=nothing
    try
        for slot in 1:workers
            put!(pool.inboxes[slot],job)
            dispatched=slot
        end
    finally
        # Always consume all acknowledgments before exposing failure or allowing
        # any callback-owned state to be reset. Each completion fits in its inbox
        # even if the controller is still submitting work to another slot.
        for slot in 1:dispatched
            result=take!(pool.completed[slot])
            isnothing(failure) && (failure=result)
        end
    end
    isnothing(failure) || throw(failure.ex)
    nothing
end

# Diagnostic timestamps share the coordinator's time origin. Each worker writes
# only its own entry, and the coordinator reads them after acknowledgment.
function dispatch_measured_merge!(job,pool,workers,t0)
    submitted=zeros(Int,workers)
    entered=similar(submitted)
    finished=similar(submitted)
    receipt_started=similar(submitted)
    received=similar(submitted)
    elapsed()=Int(time_ns()-t0)
    wrapped=slot -> begin
        entered[slot]=elapsed()
        job(slot)
        finished[slot]=elapsed()
    end
    dispatched=0
    failure=nothing
    dispatch_start=elapsed()
    dispatch_end=0
    try
        for slot in 1:workers
            submitted[slot]=elapsed()
            put!(pool.inboxes[slot],wrapped)
            dispatched=slot
        end
    finally
        dispatch_end=elapsed()
        for slot in 1:dispatched
            receipt_started[slot]=elapsed()
            result=take!(pool.completed[slot])
            received[slot]=elapsed()
            isnothing(failure) && (failure=result)
        end
    end
    isnothing(failure) || throw(failure.ex)
    Dict{String,Any}("dispatch_start_ns"=>dispatch_start,"dispatch_end_ns"=>dispatch_end,
        "submitted_ns"=>submitted,"entered_ns"=>entered,"finished_ns"=>finished,
        "receipt_started_ns"=>receipt_started,"received_ns"=>received)
end

struct ChunkScatterRun
    source_slot::Int
    source_buffer::Int
    local_range::UnitRange{Int}
    offset::Int
end

struct ChunkWave
    work::Vector{UnitRange{Int}}
    routes::Vector{Vector{ChunkScatterRun}}
end

struct ChunkThreadPlan
    workers::Int
    batches::Int
    mergers::Int
    waves::Vector{ChunkWave}
end

# Supply contiguous, approximate edge-plus-site weighted runs to individual
# workers. Spatial ownership remains based on whole cells and is independent of
# these runs. Each wave uses a bounded prefix of each worker's assigned run.
const CHUNKS_PER_WAVE=8

function chunk_thread_plan(w,workers)
    neighbours=w.neighbours
    env=neighbours.environment
    # Track active rather than padded edges without an edge-sized prefix array.
    weights=map(w.chunks) do chunk
        length(chunk.centers)+sum(chunk.centers) do ordinal
            i=neighbours.p[ordinal]
            count(view(env.active,env.offsets[1,i]:env.offsets[end,i]-1))
        end
    end
    cumulative=cumsum(weights)
    previous=0
    work=map(1:workers) do slot
        target=cld(last(cumulative)*slot,workers)
        stop=clamp(searchsortedfirst(cumulative,target),previous+1,length(cumulative)-(workers-slot))
        range=previous+1:stop
        previous=stop
        range
    end
    # Split large owner ranges even for a single group. Merge lanes can own
    # several disjoint slices when group boundaries and worker counts differ.
    target=cld(length(neighbours.p),workers)
    ranges=UnitRange{Int}[i:min(i+target-1,last(group)) for group in w.groups
                          for i in first(group):target:last(group)]
    ends=last.(ranges)
    mergers=min(length(ranges),workers)
    batches=min(CHUNKS_PER_WAVE,maximum(length,work))
    waves=ChunkWave[]
    for firstbatch in 1:batches:maximum(length,work)
        assignments=[first(r)+firstbatch-1:min(last(r),first(r)+firstbatch+batches-2) for r in work]
        routes=[ChunkScatterRun[] for _ in 1:mergers]
        # Stable source order within each wave for forces, energy and virial.
        for (slot,assignment) in enumerate(assignments), (buffer,index) in enumerate(assignment)
            chunk=w.chunks[index]
            for run in chunk.runs
                offset=chunk.destinations[first(run)]-first(run)
                lo=first(run)+offset
                hi=last(run)+offset
                part=searchsortedfirst(ends,lo)
                while lo<=hi
                    stop=min(hi,ends[part])
                    push!(routes[mod1(part,mergers)],ChunkScatterRun(slot,buffer,lo-offset:stop-offset,offset))
                    lo=stop+1
                    part+=1
                end
            end
        end
        push!(waves,ChunkWave(assignments,routes))
    end
    ChunkThreadPlan(workers,batches,mergers,waves)
end

compatible_chunk_storage(::Nothing,scatter,sums,workers,batches,capacity,limit)=false
function compatible_chunk_storage(storage,scatter,sums,workers,batches,capacity,limit)
    fields,totals,slots,wavesize,slotsize=storage[1]
    fields==scatter && totals==sums && slots>=workers && wavesize>=batches && capacity<=slotsize<=limit
end

function run_threaded_chunks!(f::F,w,requested,scatter,sums,buffers,assigned,inputs,init,timings,inline_merge) where F
    workers=min(requested,length(w.chunks))
    if isnothing(w.plan) || w.plan.workers != workers
        start=chunk_clock(timings)
        w.plan=chunk_thread_plan(w,workers)
        record_chunk_rebuild!(timings,:route_plan,start)
    end
    plan=w.plan::ChunkThreadPlan
    if isnothing(w.pool) || length(w.pool.tasks)<workers
        isnothing(w.pool) || close(w.pool)
        w.pool=ChunkWorkerPool(workers)
    end
    pool=w.pool::ChunkWorkerPool
    if !compatible_chunk_storage(pool.storage,scatter,sums,workers,plan.batches,
                                 w.observed_capacity,w.scatter_capacity)
        key=(scatter,sums,workers,plan.batches,w.observed_capacity)
        pool.storage=(key,[(map(T -> zeros(T,w.observed_capacity),scatter),map(T -> Ref(zero(T)),sums))
                           for _ in 1:plan.batches, _ in 1:workers])
    end
    scratch=[init(slot) for slot in 1:workers]
    if !isnothing(timings)
        peak=maximum(wave -> sum(length,wave.work),plan.waves)
        set_chunk_storage_stats!(timings,workers,pool.storage[1][5],peak)
    end
    run_chunk_waves!(f,w.neighbours,w.chunks,plan,pool,pool.storage[2],scratch,buffers,assigned,inputs,timings,inline_merge)
end

function run_chunk_waves!(f::F,n,chunks,plan,pool,storage,scratch,buffers,assigned,inputs,timings,inline_merge) where F
    for wave in plan.waves
        site_start=chunk_clock(timings)
        dispatch_chunk_workers!(pool,plan.workers) do slot
            start=chunk_clock(timings)
            for (buffer,index) in enumerate(wave.work[slot])
                scatter,sums=storage[buffer,slot]
                evaluate_chunk!(f,n,chunks[index],inputs,scratch[slot],scatter,sums,buffers.per_site,assigned)
            end
            record_chunk_stage!(timings,:site,slot,start)
            record_chunk_completion!(timings,slot)
        end
        record_chunk_barrier!(timings,:site,site_start)
        merge_start=chunk_clock(timings)
        trace=nothing
        if !isempty(buffers.scatter)
            merge_routes=slot -> begin
                start=chunk_clock(timings)
                for run in wave.routes[slot]
                    src=storage[run.source_buffer,run.source_slot][1]
                    merge_scatter_run!(buffers.scatter,src,run.local_range,run.offset)
                end
                record_chunk_stage!(timings,:merge,slot,start)
            end
            if inline_merge
                for slot in 1:plan.mergers
                    merge_routes(slot)
                end
            elseif isnothing(timings)
                dispatch_chunk_workers!(merge_routes,pool,plan.mergers)
            else
                trace=dispatch_measured_merge!(merge_routes,pool,plan.mergers,merge_start)
            end
        end
        sums_start=chunk_clock(timings)
        for (slot,assignment) in enumerate(wave.work), buffer in 1:length(assignment)
            merge_chunk_sums!(buffers.sum_outputs,storage[buffer,slot][2])
        end
        record_merge_wave!(timings,trace,merge_start,sums_start)
        # Every consumer has acknowledged before any of this storage is reused.
        record_chunk_slot_wait!(timings)
    end
    nothing
end

"""Optional wall-clock diagnostics for one chunk pass. Use a fresh instance per pass."""
mutable struct ChunkTimings
    site::Vector{UInt64}
    merge::Vector{UInt64}
    completed::Vector{UInt64}
    slot_wait::Vector{UInt64}
    site_barrier::UInt64
    merge_barrier::UInt64
    pass::UInt64
    export_time::UInt64
    capacity::Int
    peak_slots::Int
    merge_waves::Vector{Dict{String,Any}}
    rebuild::Dict{String,Float64}
end
ChunkTimings()=ChunkTimings(UInt64[],UInt64[],UInt64[],UInt64[],0,0,0,0,0,0,Dict{String,Any}[],Dict{String,Float64}())
@inline chunk_clock(::Nothing)=UInt64(0)
@inline chunk_clock(::ChunkTimings)=time_ns()
@inline record_chunk_stage!(::Nothing,stage,slot,start)=nothing
@inline record_chunk_stage!(t::ChunkTimings,stage,slot,start)=(getfield(t,stage)[slot]+=time_ns()-start)
@inline record_chunk_barrier!(::Nothing,stage,start)=nothing
function record_chunk_barrier!(t::ChunkTimings,stage,start)
    if stage===:site
        t.site_barrier+=time_ns()-start
    else
        t.merge_barrier+=time_ns()-start
    end
end
@inline record_chunk_completion!(::Nothing,slot)=nothing
@inline record_chunk_completion!(t::ChunkTimings,slot)=(t.completed[slot]=time_ns())
@inline record_chunk_slot_wait!(::Nothing)=nothing
function record_chunk_slot_wait!(t::ChunkTimings)
    stop=time_ns()
    for slot in eachindex(t.slot_wait); t.slot_wait[slot]+=stop-t.completed[slot]; end
end
@inline record_chunk_pass!(::Nothing,start)=nothing
@inline record_chunk_pass!(t::ChunkTimings,start)=(t.pass=time_ns()-start)
@inline record_chunk_export!(::Nothing,start)=nothing
@inline record_chunk_export!(t::ChunkTimings,start)=(t.export_time=time_ns()-start)
set_chunk_storage_stats!(::Nothing,workers,capacity,peak)=nothing
function set_chunk_storage_stats!(t::ChunkTimings,workers,capacity,peak)
    for a in (t.site,t.merge,t.completed,t.slot_wait)
        resize!(a,workers); fill!(a,0)
    end
    t.capacity=capacity
    t.peak_slots=peak
end

@inline record_merge_wave!(::Nothing,trace,start,sums_start)=nothing
function record_merge_wave!(t::ChunkTimings,trace,start,sums_start)
    stop=time_ns()
    t.merge_barrier+=stop-start
    wave=isnothing(trace) ? Dict{String,Any}() : trace
    wave["wall_ns"]=Int(stop-start)
    wave["sums_start_ns"]=Int(sums_start-start)
    wave["sums_ns"]=Int(stop-sums_start)
    push!(t.merge_waves,wave)
    nothing
end

@inline record_chunk_rebuild!(::Nothing,stage,start)=nothing
@inline record_chunk_rebuild!(t::ChunkTimings,stage,start)=(t.rebuild[string(stage)]=(time_ns()-start)*1e-9)
