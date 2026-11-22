"""Requested chunk-local scatter vectors, scalar refs and cell-ordered site views.
All fields are borrowed for the duration of a worker callback only.
"""
struct ChunkBuffers{Scatter<:NamedTuple,Sums<:NamedTuple,Sites<:NamedTuple}
    scatter::Scatter
    sum_outputs::Sums
    per_site::Sites
end

"""Owned geometry and complete directed environments, with cell-order permutation.
Rebuild after changing geometry/species. Rebuilding invalidates attached fields.
"""
struct ChunkNeighbours{E,X,H}
    environment::E
    positions::X
    cell::H
    p::Vector{Int}
    inverse::Vector{Int}
    ranges::Vector{UnitRange{Int}}
    version::UInt
end
function chunk_neighbours(sys,env,X,H;precision=eltype(X))
    w=Workspace(;width=env.width,precision,bruteforce_threshold=0)
    prepare!(w,sys,env.cutoff)
    ChunkNeighbours(env,X,H,copy(w.p),copy(w.inverse),copy(w.ranges),env.version)
end
"""prepare_neighbors(system, cutoff; width=1, precision=Float64)

Prepare fresh complete CPU neighborhoods. Species default to strict provider
identity; `channels` and `species_count` can specify an explicit model mapping.
"""
function prepare_neighbors(sys,cutoff; width=1,precision=Float64,channels=nothing,species_count=nothing)
    labels=[AtomsBase.species(sys,i) for i=1:length(sys)]
    unique_labels=unique(labels)
    ids=isnothing(channels) ? [findfirst(s -> isequal(s,l),unique_labels) for l in labels] : channels
    nz=isnothing(species_count) ? (isnothing(channels) ? length(unique_labels) : maximum(ids;init=0)) : species_count
    env=FullEnvironment(;width)
    prepare_full_environment!(env,sys,ids,nz,numeric_length(precision,cutoff);precision)
    X=reshape(precision[numeric_length(precision,AtomsBase.position(sys,i)[a]) for a=1:3,i=1:length(sys)],3,:)
    cv=AtomsBase.cell_vectors(sys); pb=AtomsBase.periodicity(sys)
    H=[pb[b] ? numeric_length(precision,cv[b][a]) : zero(precision) for a=1:3,b=1:3]
    chunk_neighbours(sys,env,X,H;precision)
end

struct SiteChunk
    centers::UnitRange{Int} # cell-order ordinals
    destinations::Vector{Int} # sorted cell-order atom indices
    runs::Vector{UnitRange{Int}} # contiguous runs in destinations (local indices)
end
mutable struct ChunkWorkgroups{N}
    neighbours::N
    scatter_capacity::Int
    chunks::Vector{SiteChunk}
    groups::Vector{UnitRange{Int}} # disjoint ranges of cell-order atoms
    fields::Dict{Symbol,Any}
    results::Any
    atom_outputs::Tuple
    guard::ReentrantLock
    active::Bool
    observed_capacity::Int
    pool::Union{Nothing,ChunkWorkerPool}
    plan::Union{Nothing,ChunkThreadPlan}
end
"""Build spatial ownership and site chunks from whole occupied cells in cell order.

`groups` requests up to that many logical owners, independently of the worker
count. `scatter_capacity` bounds each chunk footprint. An oversized single site
fails before any callback can run.
"""
function form_numa_workers(n::ChunkNeighbours,groups::Integer;scatter_capacity=4096)
    groups>0 || throw(ArgumentError("groups must be positive"))
    scatter_capacity isa Integer && scatter_capacity>0 || throw(ArgumentError("scatter_capacity must be positive"))
    chunks=SiteChunk[]
    for range in n.ranges
        firstcenter=first(range); destinations=Int[]
        for ordinal in range
            i=n.p[ordinal]
            footprint=sort!([n.inverse[j] for j in n.environment.footprints[i]])
            length(footprint)<=scatter_capacity || throw(ArgumentError("single site $i needs $(length(footprint)) scatter slots, capacity $scatter_capacity"))
            combined=union(destinations,footprint)
            if length(combined)>scatter_capacity
                push!(chunks,site_chunk(firstcenter:ordinal-1,destinations))
                firstcenter=ordinal; destinations=footprint
            else
                destinations=combined
            end
        end
        isempty(range) || push!(chunks,site_chunk(firstcenter:last(range),destinations))
    end
    owners=UnitRange{Int}[]
    cell_ends=last.(n.ranges)
    ncells=length(cell_ends)
    group_count=min(groups,ncells)
    firstcell=1
    for group in 1:group_count
        remaining=group_count-group
        target=cld(length(n.p)*group,group_count)
        lastcell=clamp(searchsortedfirst(cell_ends,target),firstcell,ncells-remaining)
        push!(owners,first(n.ranges[firstcell]):last(n.ranges[lastcell]))
        firstcell=lastcell+1
    end
    ChunkWorkgroups(n,Int(scatter_capacity),chunks,owners,Dict{Symbol,Any}(),nothing,(),
                    ReentrantLock(),false,0,nothing,nothing)
end
function site_chunk(centers,destinations)
    sort!(destinations)
    runs=UnitRange{Int}[]
    start=1
    for k=2:length(destinations)
        if destinations[k]!=destinations[k-1]+1
            push!(runs,start:k-1); start=k
        end
    end
    isempty(destinations) || push!(runs,start:length(destinations))
    SiteChunk(centers,destinations,runs)
end
function check_chunk_topology!(w)
    if w.neighbours.version!=w.neighbours.environment.version
        empty!(w.fields); w.results=nothing
        throw(ArgumentError("stale chunk topology; rebuild neighbors and workgroups"))
    end
end
struct ChunkSite{S,C,I,N}
    site::S
    chunk::C
    inputs::I
    neighbours::N
    ordinal::Int
end
for accessor in (:central_index,:central_species,:neighbour_count)
    @eval NI.$accessor(s::ChunkSite)=NI.$accessor(s.site)
end
for accessor in (:neighbour_species,:displacement,:active,:structural,:species_range)
    @eval NI.$accessor(s::ChunkSite,k)=NI.$accessor(s.site,k)
end
site_ordinal(s::ChunkSite)=s.ordinal
central_slot(s::ChunkSite)=searchsortedfirst(s.chunk.destinations,s.neighbours.inverse[NI.central_index(s)])
function neighbor_slot(s::ChunkSite,k)
    NI.active(s,k) || return 0
    j=s.site.environment.neighbours[edge_index(s.site,k)]
    searchsortedfirst(s.chunk.destinations,s.neighbours.inverse[j])
end
central_feature(s::ChunkSite,::Val{name}) where name=getproperty(s.inputs,name)[s.neighbours.inverse[NI.central_index(s)]]
function neighbor_feature(s::ChunkSite,k,::Val{name}) where name
    NI.active(s,k) || throw(ArgumentError("inactive edges have no neighbor feature"))
    j=s.site.environment.neighbours[edge_index(s.site,k)]
    getproperty(s.inputs,name)[s.neighbours.inverse[j]]
end
struct ChunkWorker{N,I,S,B}
    neighbours::N
    inputs::I
    scratch::S
    buffers::B
    chunk::SiteChunk
end
get_buffers(w::ChunkWorker)=w.buffers
get_scratch(w::ChunkWorker)=w.scratch
Base.length(w::ChunkWorker)=length(w.chunk.centers)
function Base.iterate(w::ChunkWorker,k=1)
    k>length(w) && return nothing
    i=w.neighbours.p[w.chunk.centers[k]]
    site=SiteView(w.neighbours.environment,w.neighbours.positions,w.neighbours.cell,i)
    ChunkSite(site,w.chunk,w.inputs,w.neighbours,k),k+1
end

# A site result is published only when every center has been assigned, including
# isbits values for which `isassigned` on ordinary arrays cannot detect omissions.
struct AssignedSiteVector{T,V<:AbstractVector{T},M} <: AbstractVector{T}
    values::V
    assigned::M
end
Base.size(v::AssignedSiteVector)=size(v.values)
Base.IndexStyle(::Type{<:AssignedSiteVector})=IndexLinear()
Base.getindex(v::AssignedSiteVector,i::Int)=v.values[i]
function Base.setindex!(v::AssignedSiteVector,x,i::Int)
    v.values[i]=x; v.assigned[i]=true
    x
end

function result_names(scatter,sums,sites)
    names=(keys(scatter)...,keys(sums)...,keys(sites)...)
    length(unique(names))==length(names) || throw(ArgumentError("ambiguous result names"))
    names
end
# The lock protects publication as well as pass/close operations. The active bit
# rejects reentry by the owning task, since the lock itself is reentrant.
function acquire_chunk_access!(guard,active,description)
    trylock(guard) || throw(ArgumentError(description))
    if active()
        unlock(guard)
        throw(ArgumentError(description))
    end
    nothing
end

"""Run a complete pass with exclusive scratch and bounded scatter buffers.

One worker runs directly on the caller. Multiple workers use a persistent Base
pool. `init_scratch()` supplies fresh exclusive scratch for each worker, once
per pass. With `slot_scratch=true`, the initializer receives the 1-based slot
index, which allows the caller to reuse its own scratch across passes. A callback
handles one chunk. Borrowed buffers and inputs are valid only during that call.
"""
function foreach_worker!(f::F,w::ChunkWorkgroups,workers::Integer;
                         scatter=NamedTuple(),sum_outputs=NamedTuple(),per_site=NamedTuple(),
                         site_inputs=(),init_scratch=()->nothing,slot_scratch=false,timings=nothing) where F
    workers>0 || throw(ArgumentError("workers must be positive"))
    acquire_chunk_access!(w.guard,()->w.active,"workgroup already has an active pass")
    w.active=true
    try
        check_chunk_topology!(w)
        names=result_names(scatter,sum_outputs,per_site)
        length(unique(site_inputs))==length(site_inputs) || throw(ArgumentError("duplicate site inputs"))
        isempty(intersect(site_inputs,names)) || throw(ArgumentError("cannot read and write the same site field"))
        all(k -> haskey(w.fields,k),site_inputs) || throw(ArgumentError("missing or incomplete site input"))
        inputs=NamedTuple{Tuple(site_inputs)}(Tuple(w.fields[k] for k in site_inputs))
        w.results=nothing
        w.observed_capacity=maximum((length(c.destinations) for c in w.chunks);init=0)
        n=length(w.neighbours.p)
        final_scatter=map(T -> zeros(T,n),scatter)
        totals=map(T -> Ref(zero(T)),sum_outputs)
        final_sites=map(T -> Vector{T}(undef,n),per_site)
        # Chunks write disjoint elements. A BitVector would share words at
        # chunk boundaries and could lose concurrent assignment flags.
        assigned=map(_ -> zeros(Bool,n),per_site)
        buffers=ChunkBuffers(final_scatter,totals,final_sites)
        init=slot_scratch ? init_scratch : _ -> init_scratch()
        start=chunk_clock(timings)
        if workers==1 || length(w.chunks)==1
            run_serial_chunks!(f,w,scatter,sum_outputs,buffers,assigned,inputs,init(1),timings)
        elseif !isempty(w.chunks)
            run_threaded_chunks!(f,w,Int(workers),scatter,sum_outputs,buffers,assigned,inputs,init,timings)
        end
        all(all,assigned) || throw(ArgumentError("incomplete per-site output; assign every center before publishing"))
        check_chunk_topology!(w)
        for name in keys(per_site); w.fields[name]=getproperty(final_sites,name); end
        w.atom_outputs=(keys(scatter)...,keys(per_site)...)
        w.results=merge(final_scatter,map(r -> r[],totals),final_sites)
        record_chunk_pass!(timings,start)
    finally
        w.active=false
        unlock(w.guard)
    end
    w
end

@inline function evaluate_chunk!(f,n,chunk,inputs,scratch,storage,local_sums,final_sites,assigned)
    active_scatter=map(storage) do array
        v=view(array,1:length(chunk.destinations)); fill!(v,zero(eltype(v))); v
    end
    for ref in local_sums; ref[]=zero(typeof(ref[])); end
    active_sites=map(final_sites,assigned) do a,mask
        AssignedSiteVector(view(a,chunk.centers),view(mask,chunk.centers))
    end
    f(ChunkWorker(n,inputs,scratch,ChunkBuffers(active_scatter,local_sums,active_sites),chunk))
    nothing
end

@inline function merge_scatter_run!(destinations,storage,local_range,offset)
    for name in keys(destinations)
        dest=getproperty(destinations,name); src=getproperty(storage,name)
        for k in local_range; dest[k+offset]+=src[k]; end
    end
end
@inline function merge_chunk_sums!(totals,local_sums)
    for name in keys(totals); getproperty(totals,name)[]+=getproperty(local_sums,name)[]; end
end

function run_serial_chunks!(f::F,w,scatter,sums,buffers,assigned,inputs,scratch,timings) where F
    set_chunk_storage_stats!(timings,1,w.observed_capacity,isempty(w.chunks) ? 0 : 1)
    storage=map(T -> zeros(T,w.observed_capacity),scatter)
    local_sums=map(T -> Ref(zero(T)),sums)
    for chunk in w.chunks
        start=chunk_clock(timings)
        evaluate_chunk!(f,w.neighbours,chunk,inputs,scratch,storage,local_sums,buffers.per_site,assigned)
        record_chunk_stage!(timings,:site,1,start)
        start=chunk_clock(timings)
        for run in chunk.runs
            merge_scatter_run!(buffers.scatter,storage,run,chunk.destinations[first(run)]-first(run))
        end
        merge_chunk_sums!(buffers.sum_outputs,local_sums)
        record_chunk_stage!(timings,:merge,1,start)
    end
end

"""Owning results; per-atom fields are exported in original atom order."""
function get_results(w::ChunkWorkgroups)
    acquire_chunk_access!(w.guard,()->w.active,"pass has not completed")
    try
        check_chunk_topology!(w)
        isnothing(w.results) && throw(ArgumentError("no completed pass"))
        NamedTuple{keys(w.results)}(Tuple(k in w.atom_outputs ? w.results[k][w.neighbours.inverse] : w.results[k] for k in keys(w.results)))
    finally
        unlock(w.guard)
    end
end
function export_field(w::ChunkWorkgroups,name::Symbol)
    acquire_chunk_access!(w.guard,()->w.active,"pass has not completed")
    try
        check_chunk_topology!(w)
        haskey(w.fields,name) || throw(ArgumentError("missing or incomplete site field $name"))
        w.fields[name][w.neighbours.inverse]
    finally
        unlock(w.guard)
    end
end

"""Refresh original-order positions for a graph whose topology is still valid.

The caller guarantees unchanged atom identities, cell mapping and complete
periodic image graph. Existing cell-ordered fields remain available. Rebuild
neighbors and workgroups if the graph or mapping changes.
"""
function refresh_positions!(w::ChunkWorkgroups,positions::AbstractMatrix)
    acquire_chunk_access!(w.guard,()->w.active,"pass has not completed")
    try
        check_chunk_topology!(w)
        X=w.neighbours.positions
        axes(positions)==axes(X) || throw(DimensionMismatch("chunk positions"))
        converted=eltype(positions)===eltype(X) ? positions : Matrix{eltype(X)}(positions)
        all(isfinite,converted) || throw(ArgumentError("nonfinite chunk positions"))
        copyto!(X,converted)
        w.results=nothing
    finally
        unlock(w.guard)
    end
    w
end

"""Wake and join the idle workers; the next pass can start a new pool."""
function Base.close(w::ChunkWorkgroups)
    acquire_chunk_access!(w.guard,()->w.active,"cannot close a workgroup during an active pass")
    try
        isnothing(w.pool) || close(w.pool)
        w.pool=nothing
    finally
        unlock(w.guard)
    end
    nothing
end

"""Route a local site callback output using the shared destination map and the
existing +g(center), -g(neighbor), -d⊗g virial convention."""
function scatter_site_output!(buffers,site,out,props)
    wants(AC.Energy,props) && (buffers.sum_outputs.energy[]+=NI.output_energy(out))
    (wants(AC.Forces,props) || wants(AC.Virial,props)) || return nothing
    for k=1:NI.neighbour_count(site)
        NI.active(site,k) || continue
        g=NI.output_gradient(out,k)
        if wants(AC.Forces,props)
            buffers.scatter.forces[central_slot(site)]+=g
            buffers.scatter.forces[neighbor_slot(site,k)]-=g
        end
        if wants(AC.Virial,props)
            d=NI.displacement(site,k)
            buffers.sum_outputs.virial[]-=d*transpose(g)
        end
    end
    nothing
end

# Generic direct/batch callback adoption for the chunk executor.
function evaluate_chunk_sites!(worker,prepared,ps,scratch,outputs,ordinals,centers,props,selected)
    count=0
    # Existing batch callbacks remain usable; only a bounded batch of local
    # outputs is retained, consumed before its next reuse.
    sites=ordinals
    for site in worker
        selected[NI.central_index(site)] || continue
        if length(outputs)==1
            NI.evaluate_site!(outputs[1],prepared,site,ps,scratch,props)
            scatter_site_output!(get_buffers(worker),site,outputs[1],props)
        else
            count+=1; sites[count]=site_ordinal(site)
            if count==length(outputs)
                evaluate_chunk_batch!(worker,prepared,ps,scratch,outputs,centers,props,sites,count); count=0
            end
        end
    end
    count>0 && evaluate_chunk_batch!(worker,prepared,ps,scratch,outputs,centers,props,sites,count)
    nothing
end
function evaluate_chunk_batch!(worker,prepared,ps,scratch,outputs,centers,props,ordinals,count)
    n=worker.neighbours
    for k=1:count; centers[k]=n.p[worker.chunk.centers[ordinals[k]]]; end
    batch=SiteBatch(n.environment,n.positions,n.cell,centers,count)
    NI.evaluate_sites!(outputs,prepared,batch,ps,scratch,props)
    for k=1:count
        ordinal=ordinals[k]
        site=ChunkSite(NI.site_at(batch,k),worker.chunk,worker.inputs,n,ordinal)
        scatter_site_output!(get_buffers(worker),site,outputs[k],props)
    end
    nothing
end
