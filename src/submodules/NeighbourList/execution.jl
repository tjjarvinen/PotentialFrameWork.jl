workspace_width(::Workspace{W}) where W=W

# Request normalization belongs at the boundary, never inside a model callback.
function normalize_properties(properties)
    props=properties isa Tuple ? properties : (properties,)
    all(p -> p isa Union{AC.Energy,AC.Forces,AC.Virial},props) ||
        throw(ArgumentError("supported properties are Energy(), Forces() and Virial()"))
    e=any(p -> p isa AC.Energy,props)
    f=any(p -> p isa AC.Forces,props)
    v=any(p -> p isa AC.Virial,props)
    ((e ? (AC.Energy(),) : ())...,(f ? (AC.Forces(),) : ())...,(v ? (AC.Virial(),) : ())...)
end
@inline wants(::Type{P},props::Tuple) where P = any(p -> p isa P,props)
function check_units(p)
    AC.length_unit(p)==u"Å" && AC.energy_unit(p)==u"eV" ||
        throw(ArgumentError("native CPU execution requires declared Å/eV units"))
    nothing
end

NI.central_species(b::PairBlock)=b.central_species
NI.neighbour_species(b::PairBlock)=b.neighbour_species
NI.neighbour_species(b::PairBlock,k)=b.neighbour_species[k]
NI.neighbour_count(::PairBlock{W}) where W=W
NI.displacement(b::PairBlock)=(b.dx,b.dy,b.dz)
@inline NI.distance(b::PairBlock)=sqrt(b.dx*b.dx+b.dy*b.dy+b.dz*b.dz)
NI.active(b::PairBlock)=b.active
NI.active(b::PairBlock,k)=lane(b.active,k)
NI.structural(b::PairBlock)=b.structural
NI.structural(b::PairBlock,k)=lane(b.structural,k)
mutable struct PairOutput{V}
    energy::V
    gradient::NTuple{3,V}
    radial_derivative::V
    radial::Bool
end
# Retain the existing Cartesian output constructor.
PairOutput(e,g::NTuple{3,V}) where V=PairOutput(e,g,zero(e),false)
@inline NI.set_energy!(o::PairOutput,e)=(o.energy=e; nothing)
@inline NI.set_gradient!(o::PairOutput,g)=(o.gradient=g; o.radial=false; nothing)
@inline NI.set_radial_derivative!(o::PairOutput,g)=(o.radial_derivative=g; o.radial=true; nothing)
NI.output_radial_derivative(o::PairOutput)=o.radial_derivative
@inline function pair_gradient(o::PairOutput,b)
    o.radial || return o.gradient
    r=NI.distance(b)
    # Preserve the regular-origin convention: r=0 requires dE/dr=0 or rejection
    # by the radial model. Safe division also covers inactive dummy lanes.
    coefficient=o.radial_derivative/select(r>0,r,one(r))
    d=NI.displacement(b)
    ntuple(a -> coefficient*d[a],3)
end
NI.output_energy(o::PairOutput)=o.energy
NI.output_gradient(o::PairOutput)=o.gradient

struct SiteView{E,X,H}
    environment::E
    positions::X
    cell::H
    center::Int
end
NI.central_index(s::SiteView)=s.center
NI.central_species(s::SiteView)=s.environment.species[s.center]
NI.neighbour_count(s::SiteView)=s.environment.offsets[end,s.center]-s.environment.offsets[1,s.center]
@inline edge_index(s::SiteView,k)=s.environment.offsets[1,s.center]+k-1
NI.neighbour_species(s::SiteView,k)=s.environment.species[s.environment.neighbours[edge_index(s,k)]]
NI.displacement(s::SiteView,k)=edge_displacement(s.positions,s.cell,s.environment,s.center,edge_index(s,k))
NI.active(s::SiteView,k)=s.environment.active[edge_index(s,k)]
NI.structural(s::SiteView,k)=NI.active(s,k)
NI.species_range(s::SiteView,z)=(s.environment.offsets[z,s.center]-s.environment.offsets[1,s.center]+1):(s.environment.offsets[z+1,s.center]-s.environment.offsets[1,s.center])
mutable struct SiteOutput{T}
    energy::T
    gradients::Vector{SVector{3,T}}
end
SiteOutput(capacity,::Type{T}) where T=SiteOutput(zero(T),zeros(SVector{3,T},capacity))
@inline NI.set_energy!(o::SiteOutput,e)=(o.energy=e; nothing)
@inline NI.set_gradient!(o::SiteOutput,k,g)=(o.gradients[k]=g; nothing)
NI.output_energy(o::SiteOutput)=o.energy
NI.output_gradient(o::SiteOutput,k)=o.gradients[k]
struct SiteBatch{E,X,H}
    environment::E
    positions::X
    cell::H
    centers::Vector{Int}
    count::Int
end
NI.site_count(b::SiteBatch)=b.count
NI.site_at(b::SiteBatch,k)=SiteView(b.environment,b.positions,b.cell,b.centers[k])
NI.output_at(out::Vector{<:SiteOutput},k)=out[k]
# CPU loop fallback is deliberately outside the pure protocol. A model can
# specialize on its prepared evaluator to replace this one method.
function NI.evaluate_sites!(out,p,b,ps,scratch,props)
    for k=1:NI.site_count(b)
        NI.evaluate_site!(NI.output_at(out,k),p,NI.site_at(b,k),ps,scratch,props)
    end
    nothing
end
mutable struct SiteSlot{S,T}
    scratch::S
    outputs::Vector{SiteOutput{T}}
    centers::Vector{Int}
end
function site_slots(scratch,capacity,batch_size,::Type{T}) where T
    batch_size isa Integer && batch_size>0 || throw(ArgumentError("batch_size must be positive"))
    [SiteSlot(s,[SiteOutput(capacity,T) for _=1:batch_size],zeros(Int,batch_size)) for s in scratch]
end
@inline function reduce_site!(force,inverse,energies,virials,env,i,site,out,props)
    wants(AC.Energy,props) && (energies[i]=NI.output_energy(out))
    for k=1:NI.neighbour_count(site)
        NI.active(site,k) || continue
        g=NI.output_gradient(out,k)
        j=env.neighbours[edge_index(site,k)]
        if wants(AC.Forces,props)
            for a=1:3
                force[a,inverse[i]]+=g[a]; force[a,inverse[j]]-=g[a]
            end
        end
        if wants(AC.Virial,props)
            d=NI.displacement(site,k)
            for b=1:3,a=1:3; virials[a+3(b-1),i]-=d[a]*g[b]; end
        end
    end
    nothing
end
function consume_batch!(slot,count,p,env,X,H,ps,props,force,inverse,energies,virials,batch_size)
    batch=SiteBatch(env,X,H,slot.centers,count)
    if batch_size==1
        NI.evaluate_site!(slot.outputs[1],p,NI.site_at(batch,1),ps,slot.scratch,props)
    else
        NI.evaluate_sites!(slot.outputs,p,batch,ps,slot.scratch,props)
    end
    for k=1:count
        site=NI.site_at(batch,k)
        if wants(AC.Forces,props) || wants(AC.Virial,props)
            reduce_site!(force,inverse,energies,virials,env,slot.centers[k],site,slot.outputs[k],props)
        elseif wants(AC.Energy,props)
            energies[slot.centers[k]]=NI.output_energy(slot.outputs[k])
        end
    end
    nothing
end
function site_tasks!(slot,indices,mask,start,stride,p,env,X,H,ps,props,force,inverse,energies,virials,batch_size)
    count=0
    for k=start:stride:length(indices)
        i=indices[k]; mask[i] || continue
        count+=1; slot.centers[count]=i
        if count==batch_size
            consume_batch!(slot,count,p,env,X,H,ps,props,force,inverse,energies,virials,batch_size)
            count=0
        end
    end
    count>0 && consume_batch!(slot,count,p,env,X,H,ps,props,force,inverse,energies,virials,batch_size)
    nothing
end
function execute_sites!(force,inverse,energies,virials,env,selected,mask,slots,p,X,H,ps,props,batch_size;parallel=false)
    fill!(force,0); fill!(energies,0); fill!(virials,0)
    isempty(props) && return nothing
    if parallel
        for phase in env.phases
            nslots=min(length(slots),length(phase))
            Threads.@threads :dynamic for k=1:nslots
                site_tasks!(slots[k],phase,mask,k,nslots,p,env,X,H,ps,props,force,inverse,energies,virials,batch_size)
            end
        end
    else
        site_tasks!(slots[1],selected,mask,1,1,p,env,X,H,ps,props,force,inverse,energies,virials,batch_size)
    end
    nothing
end

@inline function update_virial(vir, d, f)
    ntuple(k -> vir[k]-sum(d[mod1(k,3)]*f[cld(k,3)]),Val(9))
end

@inline function cell_task!(w::Workspace{W,T}, p, req, task, parameters=nothing, scratch=nothing) where {W,T}
    e=zero(T)
    vir=ntuple(_ -> zero(T),9)
    for i in w.ranges[task]
        central=ntuple(_ -> zero(T),3)
        for image in w.images[task]
            r=w.ranges[image.cell]
            # Same-cell half-list skips the already owned lower triangle entirely.
            lo=max(first(r),i+(image.self ? 0 : 1))
            for j=lo:W:last(r)
                js,b=geometry(w,i,j,last(r),image)
                any(b.active) || continue
                out=PairOutput(zero(b.dx),(zero(b.dx),zero(b.dx),zero(b.dx)))
                NI.evaluate_pair!(out,p,b,parameters,scratch,req)
                en=select(b.active,out.energy,zero(out.energy))
                wants(AC.Energy,req) && (e+=sum(en))
                if (wants(AC.Forces,req) || wants(AC.Virial,req))
                    gradient=pair_gradient(out,b)
                    fx=select(b.active,gradient[1],zero(gradient[1]))
                    fy=select(b.active,gradient[2],zero(gradient[2]))
                    fz=select(b.active,gradient[3],zero(gradient[3]))
                    if wants(AC.Forces,req)
                        central=(central[1]+sum(fx),central[2]+sum(fy),central[3]+sum(fz))
                        # Duplicate padding destinations are never consumed.
                        for k=1:W
                            lane(b.active,k) || continue
                            dest=js[k]
                            w.force[1,dest]-=lane(fx,k)
                            w.force[2,dest]-=lane(fy,k)
                            w.force[3,dest]-=lane(fz,k)
                        end
                    end
                    f=(fx,fy,fz); d=(b.dx,b.dy,b.dz)
                    wants(AC.Virial,req) && (vir=update_virial(vir,d,f))
                end
            end
        end
        if wants(AC.Forces,req)
            for a=1:3; w.force[a,i]+=central[a]; end
        end
    end
    w.energies[task]=e
    for a=1:9; w.virials[a,task]=vir[a]; end
    nothing
end

"""Execute already prepared numeric storage. Internal benchmarking entrypoint;
positions/species/cell must not change between prepare! and this call.
Each phase is a barrier; each cell task is serial and owns its scalar partials.
"""
function execute!(w::Workspace, bound, req; parallel=false,parameters=nothing,scratch=nothing)
    w.schedule_version==w.version || error("schedule topology version mismatch")
    NI.update_prepared!(bound,parameters)
    fill!(w.force,0); fill!(w.energies,0); fill!(w.virials,0)
    isempty(req) && return w
    if parallel
        for phase in w.phases
            Threads.@threads :dynamic for k in eachindex(phase)
                cell_task!(w,bound,req,phase[k],parameters,isnothing(scratch) ? nothing : scratch[phase[k]])
            end
        end
    else
        for task in eachindex(w.ranges); cell_task!(w,bound,req,task,parameters,isnothing(scratch) ? nothing : scratch[task]); end
    end
    w
end


NI.update_prepared!(prepared,parameters)=nothing

"""Native CPU preparation for a protocol site potential, without a required
supertype. Required species are mapped explicitly; no provider ID is a channel ID.
Returns owned geometry and a prepared executor. Topology is rebuilt, never reused
with a skin. See execute! for numeric original-order inputs and ownership.
"""
function prepare_site_execution(p,sys;precision=Float64,width=1,slots=Threads.nthreads(),batch_size=1,kwargs...)
    NI.evaluation_kind(p)===:site || throw(ArgumentError("site callback required"))
    precision in (Float32,Float64) || throw(ArgumentError("CPU precision must be Float32/64"))
    slots>0 || throw(ArgumentError("slots must be positive"))
    batch_size isa Integer && batch_size>0 || throw(ArgumentError("batch_size must be positive"))
    check_units(p)
    labels=[AtomsBase.species(sys,i) for i=1:length(sys)]
    required=NI.species_channels(p)
    channels=isnothing(required) ? unique(labels) : required
    ids=map(labels) do s
        k=findfirst(t -> isequal(s,t),channels)
        isnothing(k) && throw(ArgumentError("unsupported species $s"))
        k
    end
    prepared=NI.prepare_potential(p,labels;precision,width,kwargs...)
    env=FullEnvironment(;width)
    prepare_full_environment!(env,sys,ids,length(channels),NI.search_cutoff(p))
    capacity=maximum((env.offsets[end,i]-env.offsets[1,i] for i=1:length(sys));init=0)
    scratch=[NI.make_scratch(prepared,capacity,precision) for _=1:slots]
    porder=sortperm(env.species); inverse=invperm(porder)
    execution=SiteExecution(prepared,env,site_slots(scratch,capacity,batch_size,precision),
        zeros(precision,3,length(sys)),inverse,zeros(precision,length(sys)),zeros(precision,9,length(sys)),
        falses(length(sys)),env.version,capacity,Int(batch_size),width)
    X=[ustrip(u"Å",AtomsBase.position(sys,i)[a]) for a=1:3,i=1:length(sys)]
    cv=AtomsBase.cell_vectors(sys); pb=AtomsBase.periodicity(sys)
    H=[pb[b] ? ustrip(u"Å",cv[b][a]) : zero(precision) for a=1:3,b=1:3]
    (;execution,positions=X,cell=H)
end
struct SiteExecution{P,E,S,T}
    potential::P
    environment::E
    slots::S
    force::Matrix{T}
    inverse::Vector{Int}
    energies::Vector{T}
    virials::Matrix{T}
    mask::BitVector
    version::UInt
    capacity::Int
    batch_size::Int
    width::Int
end
function selected_centers(centers,n)
    isnothing(centers) && return 1:n
    ids=collect(Int,centers)
    all(i -> 1<=i<=n,ids) || throw(ArgumentError("invalid center index"))
    length(unique(ids))==length(ids) || throw(ArgumentError("duplicate center index"))
    ids
end
"""Execute numeric current positions at fixed prepared topology. Results are
borrowed in executor storage; site_result makes owning original-order copies.
This mutating inference API does not provide new AD rules.
"""
function execute!(w::SiteExecution,X,H,ps;properties=(AC.Energy(),AC.Forces(),AC.Virial()),centers=nothing,parallel=false)
    env=w.environment
    env.version==w.version && length(env.species)==length(w.inverse) || throw(ArgumentError("stale prepared topology; prepare again"))
    env.width==w.width || throw(ArgumentError("stale prepared layout width"))
    maximum((env.offsets[end,i]-env.offsets[1,i] for i=1:length(env.species));init=0)<=w.capacity || throw(ArgumentError("stale prepared capacity"))
    size(X)==(3,length(env.species)) && size(H)==(3,3) || throw(DimensionMismatch("geometry"))
    props=normalize_properties(properties)
    selected=selected_centers(centers,length(env.species))
    fill!(w.mask,false)
    for i in selected; w.mask[i]=true; end
    NI.update_prepared!(w.potential,ps)
    execute_sites!(w.force,w.inverse,w.energies,w.virials,env,selected,w.mask,w.slots,w.potential,X,H,ps,props,w.batch_size;parallel)
    w
end
site_result(w::SiteExecution)=(;energy=sum(w.energies),forces=w.force[:,w.inverse],virial=reshape(vec(sum(w.virials;dims=2)),3,3))
