const SiteNI=PotentialFrameWork.NeighbourInterface
const SiteA=PotentialFrameWork.ACE
function SiteNI.make_scratch(p::SiteA.ScalarSiteModel,capacity,T,ps)
    D=ForwardDiff.Dual{ACEParameterDirection,T,1}
    SiteA.SiteParameterScratch(SiteA.parameter_zeros(ps,T),SiteA.parameter_zeros(ps,D))
end
function SiteNI.make_scratch(p::SiteA.FeatureEnergySiteModel,capacity,T,ps)
    D=ForwardDiff.Dual{ACEParameterDirection,T,1}
    SiteA.SiteParameterScratch(SiteA.parameter_zeros(ps.features,T),SiteA.parameter_zeros(ps.features,D))
end
product_value(x,w,directional)=directional ? w*ForwardDiff.value(x)+ForwardDiff.partials(x)[1] : w*x
function add_site_product!(b::NamedTuple,g::NamedTuple,w,directional)
    foreach(k -> add_site_product!(getproperty(b,k),getproperty(g,k),w,directional),keys(b)); nothing
end
add_site_product!(b::Base.RefValue,g::Real,w,directional)=(b[]+=product_value(g,w,directional); nothing)
function add_site_product!(b::AbstractArray,g::AbstractArray,w,directional)
    for k in eachindex(b,g); b[k]+=product_value(g[k],w,directional); end
    nothing
end

function SiteA.feature_radial_parameters_vjp(l::SiteA.ACEFeatureLayer,data,iz,ps::NamedTuple,adj)
    for z in unique(data.zs); SiteA.feature_require_parameter_ad(l.radial[iz,z]); end
    flat=parameter_flatten(ps)
    g=ForwardDiff.gradient(flat) do values
        params=parameter_rebuild(ps,values)
        sum((dot(view(adj,j,:),SiteA.radial_values(l.radial[iz,data.zs[j]],norm(data.Rs[j]),params)) for j in eachindex(data.Rs));init=zero(eltype(values)))
    end
    parameter_rebuild(ps,g)
end

# Included in the optional ForwardDiff extension. Local head AD is independent
# of feature dimension, neighbor topology and ACE reverse implementation.
function PotentialFrameWork.ACE.head_reverse(head,h,ps;parameters=true)
    nh=sum(length,h)
    hflat=vcat(map(vec,h)...)
    flat=parameters ? vcat(hflat,parameter_flatten(ps)) : hflat
    function unpack_h(x)
        offset=0
        map(h) do block
            n=length(block); lo=offset+1; offset+=n
            reshape(x[lo:offset],size(block))
        end
    end
    g=ForwardDiff.gradient(flat) do x
        hp=parameters ? parameter_rebuild(ps,x[nh+1:end]) : ps
        head(unpack_h(x),hp)
    end
    (;energy=head(h,ps),features=unpack_h(g),
      parameters=parameters ? parameter_rebuild(ps,g[nh+1:end]) : nothing)
end

function PotentialFrameWork.ACE.feature_energy_product(p,ps,iz,zs,Rs,edges,de,dirs;
        request=PotentialFrameWork.NeighbourInterface.SiteInputs(),parameter_bar=nothing,scratch=nothing)
    A=PotentialFrameWork.ACE
    directional=!isnothing(dirs) && any(v -> any(!iszero,v),dirs)
    # Derivatives through gradients need generic second/mixed radial arithmetic,
    # including parameter-only products; first-order features alone do not.
    if directional
        for ch in p.features.layer.radial; A.feature_require_parameter_ad(ch); end
    end
    localrs=directional ? [SVector{3}(ntuple(a -> ForwardDiff.Dual{ACEParameterDirection}(Rs[j][a],dirs[j][a]),3)) for j in eachindex(Rs)] : Rs
    data=A.feature_local_data(p.features.layer,ps.features,iz,zs,localrs,edges,p.features.basis)
    buffer=isnothing(scratch) ? nothing : directional ? scratch.directional : scratch.plain
    rev=A.feature_energy_reverse(p,ps,iz,data;request,gradient=buffer)
    value(x)=directional ? ForwardDiff.value(x) : x
    tangent(x)=directional ? ForwardDiff.partials(x)[1] : zero(x)
    weight=isnothing(de) ? 0 : de
    dr=request.displacements ? [map(x -> weight*value(x)+tangent(x),g) for g in rev.displacements] : nothing
    gp=if request.parameters
        if isnothing(parameter_bar)
            g=parameter_flatten(rev.parameters)
            parameter_rebuild(ps,map(x -> weight*value(x)+tangent(x),g))
        else
            add_site_product!(parameter_bar,rev.parameters,weight,directional)
            nothing
        end
    else
        nothing
    end
    (;displacements=dr,parameters=gp)
end

function PotentialFrameWork.ACE.feature_jvp(X,l,ps,env,H,V;parameter_direction=nothing)
    size(V)==size(X) || throw(DimensionMismatch("feature position direction"))
    dualX=map((x,v) -> ForwardDiff.Dual{ACEParameterDirection}(x,v),X,V)
    dualps=if isnothing(parameter_direction)
        ps
    else
        PotentialFrameWork.ACE.check_parameter_direction(ps,parameter_direction)
        flat=parameter_flatten(ps); direction=parameter_flatten(parameter_direction)
        length(flat)==length(direction) || throw(DimensionMismatch("feature parameter direction"))
        parameter_rebuild(ps,map((x,v) -> ForwardDiff.Dual{ACEParameterDirection}(x,v),flat,direction))
    end
    h=PotentialFrameWork.ACE.ace_features(dualX,l,dualps,env,H)
    map(block -> map(x -> ForwardDiff.partials(x)[1],block),h)
end
