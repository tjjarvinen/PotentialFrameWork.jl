"""Small CPU composition boundary for equivariant features and a smooth invariant
head. `head(blocks, head_parameters)` returns one scalar site energy using generic
arithmetic. Fixed head code/configuration belongs in `head`; all trainable numbers
belong in `(features=feature_ps, head=head_ps)`. Optional ForwardDiff supplies
local head reverse/mixed products. This is not a message-passing model.
"""
struct ACEFeatureEnergy{L,H}
    layer::L
    head::H
end
struct FeatureEnergySiteModel{F,H}
    features::F
    head::H
end
FeatureEnergySiteModel(m::ACEFeatureEnergy)=FeatureEnergySiteModel(FeatureSiteModel(m.layer),m.head)
function head_reverse end
function feature_energy_product end
NI.make_scratch(::FeatureEnergySiteModel,capacity,T)=nothing
NI.site_numeric_type(::FeatureEnergySiteModel,T)=promote_type(Float64,T)
NI.make_site_output(p::FeatureEnergySiteModel,capacity,T,outputs)=NI.make_site_output(p.features,capacity,T,outputs)
NI.site_derivative_capabilities(::FeatureEnergySiteModel)=(energy=NI.SiteInputs(),gradients=NI.SiteInputs(),features=NI.SiteInputs())
function feature_readout(l,ps,iz,data)
    Tuple(data.B[b]*view(getproperty(ps.mixing,feature_key(b)),:,:,iz) for b in eachindex(l.irreps))
end
function feature_energy_reverse(p,ps,iz,data;request=NI.SiteInputs(),gradient=nothing)
    h=feature_readout(p.features.layer,ps.features,iz,data)
    head=head_reverse(p.head,h,ps.head;parameters=request.parameters)
    reverse=feature_local_reverse(p.features.layer,ps.features,iz,data,head.features;request,gradient)
    (;energy=head.energy,features=h,displacements=reverse.displacements,
      parameters=(features=reverse.parameters,head=head.parameters))
end
function NI.evaluate_site!(out,p::FeatureEnergySiteModel,site,ps,scratch,outputs::NI.SiteOutputs)
    T=eltype(out.features[1])
    data,h=feature_site_forward(p.features,site,ps.features,T)
    outputs.energy && (out.energy=p.head(h,ps.head))
    if outputs.gradients
        rev=feature_energy_reverse(p,ps,NI.central_species(site),data;request=NI.SiteInputs(true,false))
        for (j,k) in enumerate(data.edges); out.gradients[k]=rev.displacements[j]; end
    end
    if outputs.features
        for b in eachindex(h); out.features[b].=h[b]; end
    end
    nothing
end
function NI.site_vjp!(rb,pb,p::FeatureEnergySiteModel,site,ps,bars,scratch;request=NI.SiteInputs())
    T=eltype(eltype(rb)); edges,Rs,zs=local_geometry(site,T)
    dirs=isnothing(bars.gradients) ? nothing : [bars.gradients[k] for k in edges]
    if !isnothing(bars.energy) || !isnothing(bars.gradients)
        g=feature_energy_product(p,ps,NI.central_species(site),zs,Rs,edges,bars.energy,dirs;request,parameter_bar=request.parameters ? pb : nothing,scratch)
        if request.displacements
            for (j,k) in enumerate(edges); rb[k]+=g.displacements[j]; end
        end
    end
    if !isnothing(bars.features)
        data=feature_local_data(p.features.layer,ps.features,NI.central_species(site),zs,Rs,edges,p.features.basis)
        fg=feature_local_reverse(p.features.layer,ps.features,NI.central_species(site),data,bars.features;request,gradient=scratch.plain)
        if request.displacements
            for (j,k) in enumerate(edges); rb[k]+=fg.displacements[j]; end
        end
        request.parameters && NI.add_parameter_buffer!(p.features,pb.features,fg.parameters)
    end
    nothing
end
"""Evaluate a nonlinear invariant feature head using generic serial site chunks.
Returns numeric E/F/W; optional `outputs=SiteOutputs(true,true,true)` also exports
original-order model feature values. Energy/force/virial losses differentiate
into feature parameters, head parameters and positions at fixed cell/topology.
No arbitrary nested reverse AD is implied.
"""
function feature_energy_evaluate(X,m::ACEFeatureEnergy,ps,n;kwargs...)
    NL.site_evaluate(X,FeatureEnergySiteModel(m),ps,n;kwargs...)
end

"""Feature position JVP in the public (m,channel,atom) block layout. No Jacobian
is formed; positions/topology/cell have the ace_features conventions."""
function feature_jvp end
NL.prepare_neighbors(sys,m::ACEFeatureEnergy;kwargs...)=NL.prepare_neighbors(sys,m.layer;kwargs...)

function check_parameter_direction(ps::NamedTuple,v::NamedTuple)
    keys(ps)==keys(v) || throw(ArgumentError("parameter direction field names/order"))
    foreach(k -> check_parameter_direction(getproperty(ps,k),getproperty(v,k)),keys(ps))
end
check_parameter_direction(ps::AbstractArray{<:Real},v::AbstractArray{<:Real})=
    axes(ps)==axes(v) ? nothing : throw(DimensionMismatch("parameter direction axes"))
check_parameter_direction(ps::Real,v::Real)=nothing
check_parameter_direction(ps,v)=throw(ArgumentError("parameter direction leaf type/structure"))
