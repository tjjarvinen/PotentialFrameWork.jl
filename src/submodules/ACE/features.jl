"""CPU per-atom equivariant features, separate from scalar energies.

Construct with `ACEFeatureLayer(; species, radial, correlation, lmax,
output_irreps, cutoff=nothing, totaldegree=...)`. `radial` is a CallableRadial
or a central×neighbor matrix of supported radial channels; its (n,l) labels
specify the density factors. `totaldegree` bounds sum(n+l). See
`docs/ace_features.md` for conventions, resolved construction and derivatives.
"""
struct ACEFeatureLayer{R,I,C}
    atomic_numbers::Vector{Int}
    radial::R
    radial_spec::Vector{Tuple{Int,Int}}
    lmax::Int
    density::Vector{NTuple{3,Int}} # neighbor channel (0 = summed), radial, harmonic
    products::Vector{Vector{Int}}
    coupling::C                  # m-fast, basis channel, product
    irreps::I
    cutoff::Float64
end

function ACEFeatureLayer(;species,radial,correlation,lmax,output_irreps,cutoff=nothing,totaldegree=nothing,radial_spec=nothing)
    correlation isa Integer && !(correlation isa Bool) && correlation>0 || throw(ArgumentError("correlation must be a positive integer"))
    lmax isa Integer && !(lmax isa Bool) && lmax>=0 || throw(ArgumentError("lmax must be a nonnegative integer"))
    irreps=feature_irreps(output_irreps)
    z=[Int(AtomsBase.atomic_number(ChemicalSpecies(s))) for s in species]
    !isempty(z) && length(unique(z))==length(z) || throw(ArgumentError("species must be nonempty and unique"))
    all(isequal(ChemicalSpecies(s),ChemicalSpecies(zi)) for (s,zi) in zip(species,z)) || throw(ArgumentError("feature species must be elements without isotope/name labels"))
    nz=length(z)
    radial=radial isa AbstractMatrix ? deepcopy(radial) : fill(deepcopy(radial),nz,nz)
    size(radial)==(nz,nz) || throw(DimensionMismatch("radial central×neighbor channels"))
    rspec=isnothing(radial_spec) ? channel_spec(first(radial)) : Tuple{Int,Int}.(radial_spec)
    !isempty(rspec) && all(n>0 && l>=0 for (n,l) in rspec) && length(unique(rspec))==length(rspec) ||
        throw(ArgumentError("radial spec must have unique positive n, nonnegative l labels"))
    all(ch isa RadialChannel ? size(ch.coefficients,2)==length(rspec) : channel_spec(ch)==rspec for ch in radial) || throw(ArgumentError("radial channel ordering must agree across species"))
    rc=maximum(radial_cutoff,radial)
    if !isnothing(cutoff)
        cutoff isa Real && isfinite(cutoff) && cutoff>=rc || throw(ArgumentError("cutoff must cover every radial cutoff"))
        rc=Float64(cutoff)
    end
    td=isnothing(totaldegree) ? correlation*maximum(n+l for (n,l) in rspec) : totaldegree
    td isa Integer && !(td isa Bool) && td>0 || throw(ArgumentError("totaldegree must be a positive integer"))
    density,products,C=generate_feature_tables(nz,rspec,Int(correlation),Int(lmax),Int(td),irreps)
    ACEFeatureLayer(z,radial,copy(rspec),Int(lmax),density,products,C,irreps,Float64(rc))
end

"""Resolved numerical feature tables. `density` holds (neighbor species channel,
radial index, harmonic index); channel 0 sums species. Each coupling has
(2L+1)*basis_channels rows with m fastest. Supplied tables must be intertwiners
in the documented real convention; use generated construction for automatic
symmetry/linear-dependence handling.
"""
function ACEFeatureLayer(radial,density,products,coupling;species,radial_spec,lmax,output_irreps,cutoff=nothing)
    ir=feature_irreps(output_irreps)
    z=[Int(AtomsBase.atomic_number(ChemicalSpecies(s))) for s in species]
    !isempty(z) && length(unique(z))==length(z) || throw(ArgumentError("species must be nonempty and unique"))
    all(isequal(ChemicalSpecies(s),ChemicalSpecies(zi)) for (s,zi) in zip(species,z)) || throw(ArgumentError("feature species must be elements"))
    lmax isa Integer && !(lmax isa Bool) && lmax>=0 || throw(ArgumentError("lmax must be nonnegative"))
    r=radial isa AbstractMatrix ? deepcopy(radial) : fill(deepcopy(radial),length(z),length(z))
    size(r)==(length(z),length(z)) || throw(DimensionMismatch("radial species channels"))
    rs=Tuple{Int,Int}.(radial_spec)
    !isempty(rs) && all(n>0 && 0<=l<=lmax for (n,l) in rs) || throw(ArgumentError("invalid radial labels"))
    all(ch isa RadialChannel ? size(ch.coefficients,2)==length(rs) : channel_spec(ch)==rs for ch in r) || throw(ArgumentError("radial channel ordering"))
    d=NTuple{3,Int}.(density); p=[Int.(v) for v in products]
    all(0<=s<=length(z) && 1<=n<=length(rs) && 1<=m<=(lmax+1)^2 for (s,n,m) in d) || throw(ArgumentError("invalid density indices"))
    all(!isempty(v) && all(a -> 1<=a<=length(d),v) for v in p) || throw(ArgumentError("invalid products"))
    length(coupling)==length(ir) || throw(DimensionMismatch("coupling block count"))
    C=map(eachindex(ir)) do b
        c=SparseMatrixCSC{Float64,Int}(coupling[b]); dim=2ir[b].L+1
        size(c,2)==length(p) && size(c,1)>0 && size(c,1)%dim==0 || throw(DimensionMismatch("incomplete coupling multiplets"))
        all(isfinite,nonzeros(c)) || throw(ArgumentError("nonfinite coupling coefficients"))
        any(!iszero,nonzeros(c)) || throw(ArgumentError("no basis channels survive"))
        c
    end
    rc=maximum(radial_cutoff,r)
    if !isnothing(cutoff)
        cutoff isa Real && isfinite(cutoff) && cutoff>=rc || throw(ArgumentError("cutoff must cover radial cutoffs"))
        rc=Float64(cutoff)
    end
    ACEFeatureLayer(z,r,rs,Int(lmax),d,p,Tuple(C),ir,Float64(rc))
end

"""Use identical resolved scalar 0e basis coordinates (many-body term only).
Species-order projections are not supported by this adapter. Masked basis rows
remain zero. No reference or separate pair contribution enters features.
"""
function ACEFeatureLayer(m::ResolvedModel;channels=1)
    m.angular_kind===:real_spherical || throw(ArgumentError("resolved feature adapter requires spherical harmonics"))
    isempty(m.policy.species) && isnothing(m.policy.orders) || throw(ArgumentError("resolved feature adapter requires an unprojected basis"))
    irreps=feature_irreps(((L=0,parity=:even,channels=channels),))
    C=spdiagm(0=>collect(Float64,m.policy.mask))*m.coupling
    ACEFeatureLayer(copy(m.atomic_numbers),deepcopy(m.radial),copy(m.radial_spec),m.lmax,
        [(0,n,h) for (n,h) in m.density],deepcopy(m.products),(C,),irreps,maximum(radial_cutoff,m.radial))
end

"""Available coupled basis multiplicities and readout rank limits, in caller order."""
feature_basis_summary(l::ACEFeatureLayer)=map(eachindex(l.irreps)) do b
    s=l.irreps[b]; n=size(l.coupling[b],1)÷(2s.L+1)
    (;L=s.L,parity=s.parity,channels=s.channels,basis_channels=n,rank_limit=min(n,s.channels))
end
feature_key(b)=Symbol(:block_,b)
"""Owning explicit parameters. Mixing axes are (basis, output channel,
central species), shared across all magnetic components. Radial scales have axes
(radial channel, central species, neighbor species). Optional callable parameters
are stored under `custom_radial` and shared by callable species channels.
"""
function feature_parameters(l::ACEFeatureLayer;custom_radial=NamedTuple())
    nz=length(l.atomic_numbers)
    weights=map(eachindex(l.irreps)) do b
        n=size(l.coupling[b],1)÷(2l.irreps[b].L+1); width=l.irreps[b].channels
        W=zeros(n,width,nz)
        for z=1:nz,c=1:width; W[mod1(c,n),c,z]=1; end
        W
    end
    (;mixing=NamedTuple{Tuple(feature_key(b) for b in eachindex(weights))}(Tuple(weights)),
      radial=(;scale=ones(length(l.radial_spec),nz,nz)),custom_radial=deepcopy(custom_radial))
end

function prepare_environment!(env::FullEnvironment,sys,l::ACEFeatureLayer;kwargs...)
    channels=[element_channel(l,AtomsBase.species(sys,i)) for i=1:length(sys)]
    NL.prepare_full_environment!(env,sys,channels,length(l.atomic_numbers),l.cutoff;kwargs...)
end

function check_feature_inputs(X,l,ps,env,H)
    Set(keys(ps))==Set((:mixing,:radial,:custom_radial)) || throw(ArgumentError("feature parameters must contain mixing, radial and custom_radial"))
    keys(ps.radial)==(:scale,) || throw(ArgumentError("radial parameters must contain scale"))
    Set(keys(ps.mixing))==Set(feature_key(b) for b in eachindex(l.irreps)) || throw(ArgumentError("mixing block names"))
    size(X)==(3,length(env.species)) || throw(DimensionMismatch("positions must be 3×N"))
    size(H)==(3,3) || throw(DimensionMismatch("cell must be 3×3"))
    all(z -> 1<=z<=length(l.atomic_numbers),env.species) || throw(ArgumentError("environment species channels"))
    for b in eachindex(l.irreps)
        W=getproperty(ps.mixing,feature_key(b)); s=l.irreps[b]
        size(W)==(size(l.coupling[b],1)÷(2s.L+1),s.channels,length(l.atomic_numbers)) || throw(DimensionMismatch("mixing block $b"))
    end
    size(ps.radial.scale)==(length(l.radial_spec),length(l.atomic_numbers),length(l.atomic_numbers)) || throw(DimensionMismatch("radial scale shape"))
    nothing
end

feature_numeric_type(X,ps,H)=typeof(zero(eltype(X))+zero(eltype(H))+0.0+parameter_seed(ps))

function feature_site_data(X,l,ps,env,H,i,basis)
    edges=[k for k=env.offsets[1,i]:env.offsets[end,i]-1 if env.active[k]]
    Rs=[edge_displacement(X,H,env,i,k) for k in edges]
    zs=[env.species[env.neighbours[k]] for k in edges]; iz=env.species[i]
    T=feature_numeric_type(X,ps,H); nr=length(l.radial_spec)
    R=zeros(T,length(edges),nr); dR=similar(R); rawR=similar(R)
    for j in eachindex(edges)
        r=norm(Rs[j]); r>0 || throw(DomainError(r,"coincident physical atoms have undefined harmonics"))
        v,d=radial_ed(l.radial[iz,zs[j]],r,ps.custom_radial)
        for n=1:nr
            rawR[j,n]=v[n]; scale=ps.radial.scale[n,iz,zs[j]]
            R[j,n]=v[n]*scale; dR[j,n]=d[n]*scale
        end
    end
    Y,dY=isempty(edges) ? (zeros(T,0,(l.lmax+1)^2),zeros(SVector{3,T},0,(l.lmax+1)^2)) : angular_ed(basis,Rs)
    A=zeros(T,length(l.density))
    for (a,(z,n,m)) in enumerate(l.density),j in eachindex(edges)
        (z==0 || zs[j]==z) && (A[a]+=R[j,n]*Y[j,m])
    end
    P=[prod((A[a] for a in indices);init=one(T)) for indices in l.products]
    B=map(eachindex(l.irreps)) do b
        reshape(l.coupling[b]*P,2l.irreps[b].L+1,:)
    end
    (;edges,Rs,zs,R,dR,rawR,Y,dY,A,P,B)
end

"""`ace_features(X, layer, ps, environment, H)` returns an owning tuple of
(2L+1, channels, atoms) arrays. X is 3×N in Å; H has lattice vectors in columns.
Magnetic components follow SpheriCart real L2 order -L:L; atoms retain original
input order. Cell and full directed topology are fixed for differentiation.
"""
function ace_features(X::AbstractMatrix,l::ACEFeatureLayer,ps,env::FullEnvironment,H::AbstractMatrix)
    check_feature_inputs(X,l,ps,env,H)
    T=feature_numeric_type(X,ps,H)
    out=map(s -> zeros(T,2s.L+1,s.channels,size(X,2)),l.irreps)
    basis=angular_basis(l.lmax,Float64)
    for i in axes(X,2)
        data=feature_site_data(X,l,ps,env,H,i,basis)
        for b in eachindex(out)
            mul!(view(out[b],:,:,i),data.B[b],view(getproperty(ps.mixing,feature_key(b)),:,:,env.species[i]))
        end
    end
    out
end

# Reuse the scalar path's division-free polynomial reverse. Species-resolved
# densities have already been formed, so no selected-species projection is used.
function feature_product_reverse(l,A,weights)
    T=promote_type(eltype(A),eltype(weights)); dA=zeros(T,length(A)); dV=similar(dA)
    degree=maximum(length,l.products;init=0)
    product_projection!(zeros(T,length(l.products)),dA,dV,A,zero(A),l.products,weights,nothing,
        zeros(T,degree+1,degree+1),zeros(T,degree+1,degree+1))
    dA
end
feature_require_parameter_ad(ch::Union{RadialChannel,TabulatedRadial})=nothing
feature_require_parameter_ad(ch::SubsetRadial)=feature_require_parameter_ad(ch.parent)
feature_require_parameter_ad(ch::CallableRadial)=ch.metadata.ad ? nothing :
    throw(ArgumentError("custom radial feature parameter derivatives require ad=true and generic ForwardDiff-compatible values"))
feature_radial_parameter_vjp(ch,rs,ps,adj)=isempty(parameter_flatten(ps)) ? ps :
    throw(ArgumentError("load ForwardDiff for custom radial feature parameter derivatives"))

"""First-order serial CPU VJP; cotangents must be numerical blocks. No complete
position/parameter Jacobian is formed. ChainRules handles zero/thunked blocks.
"""
function feature_vjp(X,l,ps,env,H,cotangents)
    check_feature_inputs(X,l,ps,env,H)
    T=promote_type(feature_numeric_type(X,ps,H),map(eltype,cotangents)...)
    dx=zeros(T,size(X)); gp=parameter_rebuild(ps,zeros(T,length(parameter_flatten(ps))))
    basis=angular_basis(l.lmax,Float64)
    for i in axes(X,2)
        iz=env.species[i]; data=feature_site_data(X,l,ps,env,H,i,basis)
        weights=zeros(T,length(l.products))
        for b in eachindex(l.irreps)
            W=getproperty(ps.mixing,feature_key(b)); G=view(cotangents[b],:,:,i)
            view(getproperty(gp.mixing,feature_key(b)),:,:,iz) .+= transpose(data.B[b])*G
            dB=G*transpose(view(W,:,:,iz))
            weights .+= transpose(l.coupling[b])*vec(dB)
        end
        dA=feature_product_reverse(l,data.A,weights)
        radial_adj=zeros(T,size(data.R))
        for j in eachindex(data.edges)
            direction=data.Rs[j]/norm(data.Rs[j]); g=zero(SVector{3,T}); z=data.zs[j]
            for (a,(channel,n,m)) in enumerate(l.density)
                (channel==0 || channel==z) || continue
                adj=dA[a]
                g+=adj*(data.dR[j,n]*data.Y[j,m]*direction+data.R[j,n]*data.dY[j,m])
                gp.radial.scale[n,iz,z]+=adj*data.rawR[j,n]*data.Y[j,m]
                radial_adj[j,n]+=adj*ps.radial.scale[n,iz,z]*data.Y[j,m]
            end
            k=data.edges[j]; dest=env.neighbours[k]
            for a=1:3; dx[a,i]-=g[a]; dx[a,dest]+=g[a]; end
        end
        if !isempty(parameter_flatten(ps.custom_radial))
            for j in eachindex(data.edges)
                localgp=feature_radial_parameter_vjp(l.radial[iz,data.zs[j]],norm(data.Rs[j]),ps.custom_radial,view(radial_adj,j,:))
                # Preserve arbitrary nested scalar/array parameter structures.
                flat=parameter_flatten(gp.custom_radial)+parameter_flatten(localgp)
                gp=merge(gp,(;custom_radial=parameter_rebuild(ps.custom_radial,flat)))
            end
        end
    end
    (;positions=dx,parameters=gp)
end
