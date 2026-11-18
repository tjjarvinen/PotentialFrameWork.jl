# Reusable model parameter reverse storage; local geometry/basis temporaries are
# still allocating reference mathematics, independently of chunk scatter bounds.
struct SiteParameterScratch{P,D}
    plain::P
    directional::D
end
zero_parameter_arrays!(x::NamedTuple)=foreach(zero_parameter_arrays!,values(x))
zero_parameter_arrays!(x::AbstractArray)=(fill!(x,zero(eltype(x))); nothing)
zero_parameter_arrays!(x::Real)=nothing
parameter_zeros(ps,T)=parameter_rebuild(ps,zeros(T,length(parameter_flatten(ps))))

# Model-local adapters, shared by functional and chunk derivative paths.
struct ScalarSiteModel{M,B}
    model::M
    basis::B
end
ScalarSiteModel(m::ResolvedModel)=ScalarSiteModel(m,angular_basis(m,eltype(m.reference)))
NI.make_scratch(::ScalarSiteModel,capacity,T)=nothing
NI.make_site_output(::ScalarSiteModel,capacity,T,outputs)=NL.SiteOutput(capacity,T)
NI.site_derivative_capabilities(::ScalarSiteModel)=(energy=NI.SiteInputs(),gradients=NI.SiteInputs(),features=NI.SiteInputs(false,false))
function local_geometry(site,::Type{T}) where T
    edges=[k for k=1:NI.neighbour_count(site) if NI.active(site,k)]
    Rs=SVector{3,T}[NI.displacement(site,k) for k in edges]
    zs=[NI.neighbour_species(site,k) for k in edges]
    edges,Rs,zs
end
function NI.evaluate_site!(out,p::ScalarSiteModel,site,ps,scratch,outputs::NI.SiteOutputs)
    outputs.features && throw(ArgumentError("scalar ACE has no feature outputs"))
    T=typeof(zero(eltype(eltype(out.gradients)))+parameter_seed(ps))
    edges,Rs,zs=local_geometry(site,T)
    data=site_data(p.model,ps,NI.central_species(site),zs,Rs;basis=p.basis)
    outputs.energy && NI.set_energy!(out,data.energy)
    if outputs.gradients
        rev=site_reverse(p.model,ps,NI.central_species(site),data)
        for (j,k) in enumerate(edges); NI.set_gradient!(out,k,rev.gradients[j]); end
    end
    nothing
end
function scalar_site_product end # implemented by optional ForwardDiff extension
function NI.site_vjp!(rb,pb,p::ScalarSiteModel,site,ps,bars,scratch;request=NI.SiteInputs())
    T=eltype(eltype(rb)); edges,Rs,zs=local_geometry(site,T)
    dirs=isnothing(bars.gradients) ? nothing : [bars.gradients[k] for k in edges]
    g=scalar_site_product(p,ps,NI.central_species(site),zs,Rs,bars.energy,dirs;request,parameter_bar=request.parameters ? pb : nothing,scratch)
    if request.displacements
        for (j,k) in enumerate(edges); rb[k]+=g.displacements[j]; end
    end
    nothing
end

"""Numerical scalar ACE chunk path with explicit positions/parameters. Same E/F/W
convention and fixed-cell derivative contract as ace_evaluate. `neighbors` is
prepared with ACE species ordering; it never supplies frozen displacements."""
function ace_evaluate_chunks(X,model,ps,neighbors;kwargs...)
    r=NL.site_evaluate(X,ScalarSiteModel(model),ps,neighbors;kwargs...)
    (;energy=r.energy,forces=r.forces,virial=r.virial)
end

struct FeatureSiteModel{L,B}
    layer::L
    basis::B
end
FeatureSiteModel(l::ACEFeatureLayer)=FeatureSiteModel(l,angular_basis(l.lmax,Float64))
mutable struct FeatureSiteOutput{T,F}
    energy::T
    gradients::Vector{SVector{3,T}}
    features::F
end
NI.output_energy(o::FeatureSiteOutput)=o.energy
NI.output_gradient(o::FeatureSiteOutput,k)=o.gradients[k]
NI.output_features(o::FeatureSiteOutput)=o.features
NI.make_scratch(::FeatureSiteModel,capacity,T)=nothing
NI.make_scratch(::FeatureSiteModel,capacity,T,ps)=SiteParameterScratch(parameter_zeros(ps,T),nothing)
NI.site_numeric_type(::FeatureSiteModel,T)=promote_type(Float64,T)
function NI.make_site_output(p::FeatureSiteModel,capacity,T,outputs)
    FeatureSiteOutput(zero(T),zeros(SVector{3,T},capacity),map(s -> zeros(T,2s.L+1,s.channels),p.layer.irreps))
end
NI.site_derivative_capabilities(::FeatureSiteModel)=(energy=NI.SiteInputs(false,false),gradients=NI.SiteInputs(false,false),features=NI.SiteInputs())
function feature_site_forward(p,site,ps,T)
    edges,Rs,zs=local_geometry(site,T); iz=NI.central_species(site)
    data=feature_local_data(p.layer,ps,iz,zs,Rs,edges,p.basis)
    h=map(eachindex(p.layer.irreps)) do b
        data.B[b]*view(getproperty(ps.mixing,feature_key(b)),:,:,iz)
    end |> Tuple
    data,h
end
function NI.evaluate_site!(out,p::FeatureSiteModel,site,ps,scratch,outputs::NI.SiteOutputs)
    (outputs.energy || outputs.gradients) && throw(ArgumentError("feature-only ACE supplies no energy/gradient output"))
    data,h=feature_site_forward(p,site,ps,eltype(out.features[1]))
    for b in eachindex(h); out.features[b].=h[b]; end
    nothing
end
function NI.site_vjp!(rb,pb,p::FeatureSiteModel,site,ps,bars,scratch;request=NI.SiteInputs())
    isnothing(bars.features) && return nothing
    data,h=feature_site_forward(p,site,ps,eltype(eltype(rb)))
    g=feature_local_reverse(p.layer,ps,NI.central_species(site),data,bars.features;request,gradient=scratch.plain)
    if request.displacements
        for (j,k) in enumerate(data.edges); rb[k]+=g.displacements[j]; end
    end
    request.parameters && NI.add_parameter_buffer!(p,pb,g.parameters)
    nothing
end
function stack_feature_sites(h,l,N)
    map(eachindex(l.irreps)) do b
        s=l.irreps[b]
        isempty(h) && return zeros(2s.L+1,s.channels,N)
        cat(map(v -> v[b],h)...;dims=3)
    end |> Tuple
end
"""Chunk feature training path. Public axes match ace_features: (m,channel,atom)."""
function ace_features_chunks(X,l,ps,n;kwargs...)
    y=NL.site_evaluate(X,FeatureSiteModel(l),ps,n;outputs=NI.SiteOutputs(false,false,true),kwargs...)
    stack_feature_sites(y.features,l,size(X,2))
end

"""Prepare complete CPU neighborhoods with the feature layer's strict element
channel order and Float64 feature arithmetic."""
function NL.prepare_neighbors(sys,l::ACEFeatureLayer;precision=Float64,width=4,kwargs...)
    channels=[element_channel(l,AtomsBase.species(sys,i)) for i=1:length(sys)]
    NL.prepare_neighbors(sys,l.cutoff;precision,width,channels,species_count=length(l.atomic_numbers),kwargs...)
end
