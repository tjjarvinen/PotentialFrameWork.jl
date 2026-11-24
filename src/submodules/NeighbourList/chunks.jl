"""Requested chunk-local scatter vectors, scalar refs and cell-ordered site views.
All fields are borrowed for the duration of a worker callback only.
"""
struct ChunkBuffers{Scatter<:NamedTuple,Sums<:NamedTuple,Sites<:NamedTuple}
    scatter::Scatter
    sum_outputs::Sums
    per_site::Sites
end

# Separate build scratch from the published geometry, so a query snapshot
# never borrows either the cell workspace or the directed-edge staging arrays.
struct ChunkNeighbourBuild{W,T}
    cells::W
    cell_scratch::CellBuildScratch{T}
    coordinates::Vector{SVector{3,T}}
    winding::Vector{SVector{3,Int}}
    indices::Vector{Vector{Int32}}
    shifts::Vector{Vector{SVector{3,Int32}}}
    marks::Vector{Int}
end
function ChunkNeighbourBuild(width,::Type{T}) where T
    cells=Workspace(;width,precision=T,bruteforce_threshold=0)
    ChunkNeighbourBuild(cells,CellBuildScratch(T),SVector{3,T}[],SVector{3,Int}[],
                        Vector{Int32}[],Vector{SVector{3,Int32}}[],Int[])
end

"""Owned geometry and complete directed environments, with cell-order permutation."""
mutable struct ChunkNeighbours{E,X,H,B}
    environment::E
    positions::X
    cell::H
    p::Vector{Int}
    inverse::Vector{Int}
    ranges::Vector{UnitRange{Int}}
    version::UInt
    build::B
    valid::Bool
    writing::Bool
    readers::Int
    guard::ReentrantLock
end
function ChunkNeighbours(env,X,H,p,inverse,ranges,version)
    build=ChunkNeighbourBuild(env.width,eltype(X))
    ChunkNeighbours(env,X,H,p,inverse,ranges,version,build,true,false,0,ReentrantLock())
end

# Compatibility helper for callers that already own a packed environment.
function chunk_neighbours(sys,env,X,H;precision=eltype(X))
    w=Workspace(;width=env.width,precision,bruteforce_threshold=0)
    prepare!(w,sys,env.cutoff)
    ChunkNeighbours(env,X,H,w.p,w.inverse,w.ranges,env.version)
end

"""Prepare complete CPU neighborhoods with owned geometry and rebuild storage.

Species default to strict provider identity. Supply `channels` and
`species_count` for model-specific ordering. Set `coloring=false` for chunk
consumers that only need footprints and cell order.
"""
function prepare_neighbors(sys,cutoff; width=1,precision=Float64,kwargs...)
    env=FullEnvironment(;width)
    build=ChunkNeighbourBuild(width,precision)
    n=ChunkNeighbours(env,Matrix{precision}(undef,3,0),zeros(precision,3,3),
                      Int[],Int[],UnitRange{Int}[],env.version,build,false,false,0,ReentrantLock())
    prepare_neighbors!(n,sys,cutoff;kwargs...)
end

"""Refill `prepare_neighbors` storage with a fresh directed graph.

The neighbor object's precision and SIMD width are fixed at construction. Any
attached workgroup must be refilled after this build. Views are borrowed and
expire on rebuild. A failed build leaves this object invalid until rebuilt.
"""
function prepare_neighbors!(n::ChunkNeighbours,sys,cutoff;channels=nothing,species_count=nothing,coloring=true)
    begin_neighbor_write!(n)
    try
        n.valid=false
        env=n.environment
        env.version+=UInt(1)
        empty!(env.timings)
        refill_chunk_neighbors!(n,n.build,sys,cutoff,channels,species_count,coloring)
        n.version=env.version
        n.valid=true
    finally
        end_neighbor_write!(n)
    end
    n
end

function refill_chunk_neighbors!(n,build::ChunkNeighbourBuild{W,T},sys,cutoff,channels,species_count,coloring) where {W,T}
    env=n.environment
    cells=build.cells
    start=time_ns()
    if isnothing(channels)
        labels=[AtomsBase.species(sys,i) for i in 1:length(sys)]
        unique_labels=unique(labels)
        channels=[findfirst(s -> isequal(s,l),unique_labels) for l in labels]
    end
    nz=isnothing(species_count) ? maximum(channels;init=0) : species_count
    set_environment_species!(env,channels,length(sys),nz,numeric_length(Float64,cutoff))
    env.timings[:channels]=(time_ns()-start)*1e-9

    start=time_ns()
    cv=AtomsBase.cell_vectors(sys)
    pb=Tuple(Bool.(AtomsBase.periodicity(sys)))
    H=SMatrix{3,3,T}(ntuple(k -> length_value(T,cv[cld(k,3)][mod1(k,3)]),9))
    size(n.positions,2)==length(sys) || (n.positions=Matrix{T}(undef,3,length(sys)))
    resize!(build.coordinates,length(sys))
    for i in eachindex(build.coordinates)
        r=AtomsBase.position(sys,i)
        position=SVector{3,T}(ntuple(a -> length_value(T,r[a]),3))
        build.coordinates[i]=position
        for a in 1:3; n.positions[a,i]=position[a]; end
    end
    for b in 1:3, a in 1:3; n.cell[a,b]=pb[b] ? H[a,b] : zero(T); end
    geometry=CoordinateSystem(build.coordinates,H,pb)
    env.timings[:boundary_extraction]=(time_ns()-start)*1e-9

    # Both traversal and chunk order use this single cell generation.
    prepare!(cells,geometry,T(env.cutoff);build=build.cell_scratch,query_only=true)
    merge!(env.timings,cells.timings)
    start=time_ns()
    indices=refill_vectors!(build.indices,length(sys))
    shifts=refill_vectors!(build.shifts,length(sys))
    image_windings!(build.winding,cells,geometry)
    # Match the snapshot query's existing checked image-label range.
    length(sys)<=typemax(Int32) || throw(ArgumentError("int_type cannot represent all atom indices"))
    visit_directed_images!(cells,build.winding,Int32) do a,b,_,shift
        push!(indices[a],b); push!(shifts[a],shift)
    end
    env.timings[:directed_traversal]=(time_ns()-start)*1e-9
    pack_environment!(env,indices,shifts,nz,build.marks)
    empty!(env.phases)
    coloring && color_environment!(env)
    start=time_ns()
    for (dest,src) in ((n.p,cells.p),(n.inverse,cells.inverse),(n.ranges,cells.ranges))
        copyto!(resize!(dest,length(src)),src)
    end
    env.timings[:cell_permutation]=(time_ns()-start)*1e-9
    nothing
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
    marks::Vector{Int}
    destinations::Vector{Int}
    version::UInt
    valid::Bool
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
    w=ChunkWorkgroups(n,0,SiteChunk[],UnitRange{Int}[],Dict{Symbol,Any}(),nothing,(),
                      Int[],Int[],n.version,false,ReentrantLock(),false,0,nothing,nothing)
    form_numa_workers!(w,n,groups;scatter_capacity)
end

"""Refill a workgroup attached to `n`, keeping its idle worker pool and layout storage.

A different neighbor object requires an allocating fallback. Use the returned
workgroup and close the old one when retiring it. A refill invalidates the
previous plan, results and retained site fields. Failure leaves it invalid.
"""
function form_numa_workers!(w::ChunkWorkgroups,n::ChunkNeighbours,groups::Integer;scatter_capacity=w.scatter_capacity)
    n===w.neighbours || return form_numa_workers(n,groups;scatter_capacity)
    acquire_workgroup!(w)
    begin_neighbor_access!(w,begin_neighbor_write!)
    try
        w.valid=false
        empty!(w.fields)
        w.results=nothing
        w.atom_outputs=()
        w.plan=nothing
        check_neighbor_topology!(n)
        groups>0 || throw(ArgumentError("groups must be positive"))
        scatter_capacity isa Integer && scatter_capacity>0 || throw(ArgumentError("scatter_capacity must be positive"))
        w.scatter_capacity=Int(scatter_capacity)
        refill_site_chunks!(w,n)
        refill_owners!(w.groups,n,groups)
        w.version=n.version
        w.valid=true
    finally
        end_neighbor_write!(n)
        unlock(w.guard)
    end
    w
end

function append_site_destinations!(destinations,marks,n,footprint,stamp)
    for j in footprint
        dest=n.inverse[j]
        if marks[dest]!=stamp
            marks[dest]=stamp
            push!(destinations,dest)
        end
    end
end

function refill_site_chunks!(w,n)
    marks=resize!(w.marks,length(n.p)); fill!(marks,0)
    destinations=w.destinations
    index=0
    for range in n.ranges
        firstcenter=first(range)
        empty!(destinations)
        for ordinal in range
            i=n.p[ordinal]
            footprint=n.environment.footprints[i]
            length(footprint)<=w.scatter_capacity ||
                throw(ArgumentError("single site $i needs $(length(footprint)) scatter slots, capacity $(w.scatter_capacity)"))
            previous=length(destinations)
            append_site_destinations!(destinations,marks,n,footprint,firstcenter)
            if length(destinations)>w.scatter_capacity
                resize!(destinations,previous)
                index+=1
                refill_site_chunk!(w.chunks,index,firstcenter:ordinal-1,destinations)
                firstcenter=ordinal
                empty!(destinations)
                append_site_destinations!(destinations,marks,n,footprint,firstcenter)
            end
        end
        index+=1
        refill_site_chunk!(w.chunks,index,firstcenter:last(range),destinations)
    end
    resize!(w.chunks,index)
    nothing
end

function refill_site_chunk!(chunks,index,centers,destinations)
    if index>length(chunks)
        push!(chunks,SiteChunk(centers,Int[],UnitRange{Int}[]))
    end
    old=chunks[index]
    ds=copyto!(resize!(old.destinations,length(destinations)),destinations)
    sort!(ds)
    runs=empty!(old.runs)
    start=1
    for k=2:length(ds)
        if ds[k]!=ds[k-1]+1
            push!(runs,start:k-1); start=k
        end
    end
    isempty(ds) || push!(runs,start:length(ds))
    old.centers==centers || (chunks[index]=SiteChunk(centers,ds,runs))
    nothing
end

function refill_owners!(owners,n,groups)
    empty!(owners)
    ncells=length(n.ranges)
    group_count=min(groups,ncells)
    firstcell=1
    for group in 1:group_count
        remaining=group_count-group
        target=cld(length(n.p)*group,group_count)
        lastcell=firstcell
        while lastcell<ncells-remaining && last(n.ranges[lastcell])<target
            lastcell+=1
        end
        push!(owners,first(n.ranges[firstcell]):last(n.ranges[lastcell]))
        firstcell=lastcell+1
    end
    nothing
end

function check_neighbor_topology!(n)
    n.valid && n.version==n.environment.version ||
        throw(ArgumentError("stale chunk topology; rebuild neighbors and workgroups"))
    nothing
end

function check_chunk_topology!(w)
    n=w.neighbours
    if !w.valid || w.version!=n.version || !n.valid || n.version!=n.environment.version
        w.valid=false
        empty!(w.fields)
        w.results=nothing
        w.plan=nothing
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

# The short guard protects reader/writer registration. Independent workgroups
# can evaluate the same generation, and a writer fails while any view is borrowed.
function begin_neighbor_write!(n)
    lock(n.guard) do
        (n.writing || n.readers>0) &&
            throw(ArgumentError("cannot mutate neighbors during an active pass or rebuild"))
        n.writing=true
    end
end
function end_neighbor_write!(n)
    lock(n.guard) do
        n.writing=false
    end
end
function begin_neighbor_read!(n)
    lock(n.guard) do
        n.writing && throw(ArgumentError("neighbors are being rebuilt"))
        n.readers+=1
    end
end
function end_neighbor_read!(n)
    lock(n.guard) do
        n.readers-=1
    end
    nothing
end
function acquire_workgroup!(w)
    acquire_chunk_access!(w.guard,()->w.active,"workgroup already has an active pass")
end
function begin_neighbor_access!(w,begin_access!)
    try
        begin_access!(w.neighbours)
    catch
        unlock(w.guard)
        rethrow()
    end
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
                         site_inputs=(),init_scratch=()->nothing,slot_scratch=false,timings=nothing,_inline_merge=false) where F
    workers>0 || throw(ArgumentError("workers must be positive"))
    acquire_workgroup!(w)
    begin_neighbor_access!(w,begin_neighbor_read!)
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
            run_threaded_chunks!(f,w,Int(workers),scatter,sum_outputs,buffers,assigned,inputs,init,timings,_inline_merge)
        end
        all(all,assigned) || throw(ArgumentError("incomplete per-site output; assign every center before publishing"))
        check_chunk_topology!(w)
        for name in keys(per_site); w.fields[name]=getproperty(final_sites,name); end
        w.atom_outputs=(keys(scatter)...,keys(per_site)...)
        w.results=merge(final_scatter,map(r -> r[],totals),final_sites)
        record_chunk_pass!(timings,start)
    finally
        w.active=false
        end_neighbor_read!(w.neighbours)
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
    acquire_workgroup!(w)
    begin_neighbor_access!(w,begin_neighbor_read!)
    try
        check_chunk_topology!(w)
        isnothing(w.results) && throw(ArgumentError("no completed pass"))
        NamedTuple{keys(w.results)}(Tuple(k in w.atom_outputs ? w.results[k][w.neighbours.inverse] : w.results[k] for k in keys(w.results)))
    finally
        end_neighbor_read!(w.neighbours)
        unlock(w.guard)
    end
end
function export_field(w::ChunkWorkgroups,name::Symbol)
    acquire_workgroup!(w)
    begin_neighbor_access!(w,begin_neighbor_read!)
    try
        check_chunk_topology!(w)
        haskey(w.fields,name) || throw(ArgumentError("missing or incomplete site field $name"))
        w.fields[name][w.neighbours.inverse]
    finally
        end_neighbor_read!(w.neighbours)
        unlock(w.guard)
    end
end

"""Refresh original-order positions for a graph whose topology is still valid.

The caller guarantees unchanged atom identities, cell mapping and complete
periodic image graph. Existing cell-ordered fields remain available. Rebuild
neighbors and workgroups if the graph or mapping changes.
"""
function refresh_positions!(w::ChunkWorkgroups,positions::AbstractMatrix)
    acquire_workgroup!(w)
    begin_neighbor_access!(w,begin_neighbor_write!)
    try
        check_chunk_topology!(w)
        X=w.neighbours.positions
        axes(positions)==axes(X) || throw(DimensionMismatch("chunk positions"))
        converted=eltype(positions)===eltype(X) ? positions : Matrix{eltype(X)}(positions)
        all(isfinite,converted) || throw(ArgumentError("nonfinite chunk positions"))
        copyto!(X,converted)
        w.results=nothing
    finally
        end_neighbor_write!(w.neighbours)
        unlock(w.guard)
    end
    w
end

"""Wake and join the idle workers; the next pass can start a new pool."""
function Base.close(w::ChunkWorkgroups)
    acquire_workgroup!(w)
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
