NI.make_scratch(p,capacity,T,ps)=NI.make_scratch(p,capacity,T)
NI.site_numeric_type(p,T)=T
NI.output_features(::SiteOutput)=()
# Serial numerical derivative boundary. No model mathematics lives here.
NI.parameter_buffer(p,ps,::Type{T}) where T=parameter_buffer_tree(ps,T)
parameter_buffer_tree(ps::NamedTuple,T)=map(x -> parameter_buffer_tree(x,T),ps)
parameter_buffer_tree(x::Real,T)=Ref(zero(T))
parameter_buffer_tree(x::AbstractArray{<:Real},T)=zeros(T,size(x))
NI.reset_parameter_buffer!(p,b::NamedTuple)=(foreach(x -> NI.reset_parameter_buffer!(p,x),values(b)); nothing)
NI.reset_parameter_buffer!(p,b::Base.RefValue)=(b[]=zero(b[]); nothing)
NI.reset_parameter_buffer!(p,b::AbstractArray)=(fill!(b,zero(eltype(b))); nothing)
function NI.add_parameter_buffer!(p,b::NamedTuple,g::NamedTuple)
    keys(b)==keys(g) || throw(ArgumentError("parameter buffer field mismatch"))
    foreach(k -> NI.add_parameter_buffer!(p,getproperty(b,k),getproperty(g,k)),keys(b)); nothing
end
NI.add_parameter_buffer!(p,b::Base.RefValue,g::Real)=(b[]+=g; nothing)
NI.add_parameter_buffer!(p,b::AbstractArray,g::AbstractArray)=(axes(b)==axes(g) || throw(DimensionMismatch("parameter axes")); b .+= g; nothing)
NI.parameter_value(p,b::NamedTuple)=map(x -> NI.parameter_value(p,x),b)
NI.parameter_value(p,b::Base.RefValue)=b[]
NI.parameter_value(p,b::AbstractArray)=copy(b)

function NI.sites_vjp!(rb,pb,p,batch,ps,bars,scratch;request=NI.SiteInputs())
    for k=1:NI.site_count(batch)
        NI.site_vjp!(rb[k],pb,p,NI.site_at(batch,k),ps,bars[k],scratch;request)
    end
    nothing
end

nonzero_bar(::Nothing)=false
nonzero_bar(x::Number)=!iszero(x)
nonzero_bar(x)=any(nonzero_bar,x)
function validate_site_derivatives(p,outputs,request)
    caps=NI.site_derivative_capabilities(p)
    for key in (:energy,:gradients,:features)
        getproperty(outputs,key) || continue
        cap=getproperty(caps,key)
        for input in (:displacements,:parameters)
            getproperty(request,input) && !getproperty(cap,input) &&
                throw(ArgumentError("$(typeof(p)) does not support $key -> $input site derivatives"))
        end
    end
    nothing
end

numeric_parameter_type(x::Real)=typeof(x)
numeric_parameter_type(x::AbstractArray{<:Real})=eltype(x)
numeric_parameter_type(x::NamedTuple)=promote_type(Float32,map(numeric_parameter_type,values(x))...)
cotangent_numeric_type(::Nothing)=Float32
cotangent_numeric_type(x::Real)=typeof(x)
cotangent_numeric_type(x::AbstractArray{<:Real})=eltype(x)
cotangent_numeric_type(x)=foldl((T,v) -> promote_type(T,cotangent_numeric_type(v)),x;init=Float32)
function numerical_workgroups(X,n;scatter_capacity=4096)
    size(X)==size(n.positions) || throw(DimensionMismatch("positions must match prepared neighbors"))
    n.version==n.environment.version || throw(ArgumentError("stale topology"))
    current=ChunkNeighbours(n.environment,X,n.cell,n.p,n.inverse,n.ranges,n.version)
    form_numa_workers(current,1;scatter_capacity)
end

"""site_evaluate(X, prepared, ps, neighbors; outputs=SiteOutputs(), centers=nothing,
scatter_capacity=4096)

Owning numerical energy/forces/virial plus original-order per-atom model feature
values. Geometry is recomputed from X and saved image shifts. Fixed topology,
cell, species and model metadata are nondifferentiable. Serial chunks; no mutable
executor is a tape. A custom ChainRules boundary owns only mathematical inputs.
"""
function site_evaluate(X,p,ps,n::ChunkNeighbours;outputs=NI.SiteOutputs(),centers=nothing,scatter_capacity=4096,batch_size=1,input_derivatives=NI.SiteInputs())
    batch_size isa Integer && batch_size>0 || throw(ArgumentError("batch_size must be positive"))
    T=NI.site_numeric_type(p,promote_type(eltype(X),numeric_parameter_type(ps)))
    w=numerical_workgroups(X,n;scatter_capacity)
    selected=falses(size(X,2)); selected[selected_centers(centers,size(X,2))].=true
    capacity=maximum((n.environment.offsets[end,i]-n.environment.offsets[1,i] for i in axes(X,2));init=0)
    scratch=NI.make_scratch(p,capacity,T)
    out=NI.make_site_output(p,capacity,T,outputs)
    props=((outputs.energy ? (AC.Energy(),) : ())...,(outputs.gradients ? (AC.Forces(),AC.Virial()) : ())...)
    sites=outputs.features ? (features=typeof(NI.output_features(out)),) : NamedTuple()
    sums=merge(outputs.energy ? (energy=T,) : NamedTuple(),outputs.gradients ? (virial=SMatrix{3,3,T,9},) : NamedTuple())
    outs=[NI.make_site_output(p,capacity,T,outputs) for _=1:batch_size]
    foreach_worker!(w,1;scatter=outputs.gradients ? (forces=SVector{3,T},) : NamedTuple(),sum_outputs=sums,per_site=sites) do worker
        pending=Vector{typeof(first(worker))}()
        function consume!()
            batch=ChunkSiteBatch(pending)
            if batch_size==1
                NI.evaluate_site!(outs[1],p,pending[1],ps,scratch,outputs)
            else
                NI.evaluate_sites!(outs,p,batch,ps,scratch,outputs)
            end
            for (k,site) in enumerate(pending)
                scatter_site_output!(get_buffers(worker),site,outs[k],props)
                outputs.features && (get_buffers(worker).per_site.features[site_ordinal(site)]=deepcopy(NI.output_features(outs[k])))
            end
            empty!(pending)
        end
        for site in worker
            i=NI.central_index(site)
            if selected[i]
                push!(pending,site)
                length(pending)==batch_size && consume!()
            elseif outputs.features
                h=deepcopy(NI.output_features(out))
                foreach(x -> fill!(x,zero(eltype(x))),h)
                get_buffers(worker).per_site.features[site_ordinal(site)]=h
            end
        end
        isempty(pending) || consume!()
    end
    r=get_results(w)
    f=get(r,:forces,zeros(SVector{3,T},size(X,2)))
    (;energy=get(r,:energy,zero(T)),forces=reshape(collect(Iterators.flatten(f)),3,:),
      virial=get(r,:virial,zero(SMatrix{3,3,T,9})),features=get(r,:features,nothing))
end

"""Numerical chunk VJP. bars is (energy,forces,virial,features), with `nothing`
for absent fields. Additive local buffers reset once per site/pass as appropriate.
Requests can omit position Hessian products entirely. No partial result escapes.
"""
function site_pullback(X,p,ps,n,bars;outputs=NI.SiteOutputs(),centers=nothing,scatter_capacity=4096,request=NI.SiteInputs(),batch_size=1)
    batch_size isa Integer && batch_size>0 || throw(ArgumentError("batch_size must be positive"))
    size(X)==size(n.positions) || throw(DimensionMismatch("positions must match prepared neighbors"))
    n.version==n.environment.version || throw(ArgumentError("stale topology"))
    isnothing(bars.energy) || bars.energy isa Real || throw(ArgumentError("energy cotangent must be real"))
    isnothing(bars.forces) || size(bars.forces)==size(X) || throw(DimensionMismatch("force cotangent"))
    isnothing(bars.virial) || size(bars.virial)==(3,3) || throw(DimensionMismatch("virial cotangent"))
    isnothing(bars.features) || length(bars.features)==size(X,2) || throw(DimensionMismatch("feature atom cotangents"))
    selected=falses(size(X,2)); selected[selected_centers(centers,size(X,2))].=true
    active=NI.SiteOutputs(nonzero_bar(bars.energy),nonzero_bar(bars.forces)||nonzero_bar(bars.virial),nonzero_bar(bars.features))
    (!outputs.energy && active.energy || !outputs.gradients && active.gradients || !outputs.features && active.features) &&
        throw(ArgumentError("cotangent supplied for an unrequested output"))
    validate_site_derivatives(p,active,request)
    T=NI.site_numeric_type(p,promote_type(eltype(X),numeric_parameter_type(ps),cotangent_numeric_type(bars)))
    if !(active.energy || active.gradients || active.features)
        pb=NI.parameter_buffer(p,ps,T)
        NI.reset_parameter_buffer!(p,pb)
        return (;positions=request.displacements ? zeros(T,size(X)) : nothing,
                 parameters=request.parameters ? NI.parameter_value(p,pb) : nothing)
    end
    w=numerical_workgroups(X,n;scatter_capacity)
    capacity=maximum((n.environment.offsets[end,i]-n.environment.offsets[1,i] for i in axes(X,2));init=0)
    scratch=NI.make_scratch(p,capacity,T,ps)
    out=NI.make_site_output(p,capacity,T,outputs)
    pb=NI.parameter_buffer(p,ps,T)
    NI.reset_parameter_buffer!(p,pb)
    rb=zeros(SVector{3,T},capacity); gb=similar(rb)
    # Permute atom cotangents once, keeping aliased image destinations distinct.
    fb=isnothing(bars.forces) ? nothing : bars.forces[:,n.p]
    hb=isnothing(bars.features) ? nothing : bars.features[n.p]
    scatter=request.displacements ? (positions=SVector{3,T},) : NamedTuple()
    foreach_worker!(w,1;scatter) do worker
        for site in worker
            i=NI.central_index(site); selected[i] || continue
            fill!(rb,zero(eltype(rb))); fill!(gb,zero(eltype(gb)))
            if active.gradients
                for k=1:NI.neighbour_count(site)
                    NI.active(site,k) || continue
                    j=site.site.environment.neighbours[edge_index(site.site,k)]
                    b=isnothing(fb) ? zero(SVector{3,T}) : SVector{3,T}(fb[:,n.inverse[i]]-fb[:,n.inverse[j]])
                    gb[k]=isnothing(bars.virial) ? b : b-transpose(bars.virial)*NI.displacement(site,k)
                end
            end
            localbars=NI.SiteCotangents(active.energy ? bars.energy : nothing,
                active.gradients ? gb : nothing, active.features ? hb[n.inverse[i]] : nothing)
            NI.site_vjp!(rb,pb,p,site,ps,localbars,scratch;request)
            if request.displacements
                if nonzero_bar(bars.virial)
                    NI.evaluate_site!(out,p,site,ps,scratch,NI.SiteOutputs(false,true,false))
                end
                for k=1:NI.neighbour_count(site)
                    NI.active(site,k) || continue
                    g=rb[k]
                    nonzero_bar(bars.virial) && (g-=bars.virial*NI.output_gradient(out,k))
                    get_buffers(worker).scatter.positions[central_slot(site)]-=g
                    get_buffers(worker).scatter.positions[neighbor_slot(site,k)]+=g
                end
            end
        end
    end
    r=get_results(w)
    dx=request.displacements ? reshape(collect(Iterators.flatten(r.positions)),3,:) : nothing
    (;positions=dx,parameters=request.parameters ? NI.parameter_value(p,pb) : nothing)
end

struct ChunkSiteBatch{V}
    sites::V
end
NI.site_count(b::ChunkSiteBatch)=length(b.sites)
NI.site_at(b::ChunkSiteBatch,k)=b.sites[k]
NI.output_at(v::AbstractVector,k)=v[k]
