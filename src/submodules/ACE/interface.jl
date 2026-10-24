# Immutable policies describe mathematics, not scratch or learned values.
struct BasisPolicy{M<:Tuple,S<:Tuple}
    origin::String
    rows::Tuple
    products::Tuple
    density::Tuple
    radial::Tuple
    mask::M
    species::S
    orders::Union{Nothing,UnitRange{Int}}
    pair::Bool
    reference::Bool
end
function default_policy(sha,rspec,density,products,C)
    id=bytes2hex(SHA.sha256(string(sha,rspec,density,products,C)))
    BasisPolicy(id,Tuple(1:size(C,1)),Tuple(eachindex(products)),Tuple(eachindex(density)),
        Tuple(eachindex(rspec)),Tuple(trues(size(C,1))),(),nothing,true,true)
end
function policy_with(p;origin=p.origin,rows=p.rows,products=p.products,density=p.density,
                     radial=p.radial,mask=p.mask,species=p.species,orders=p.orders,
                     pair=p.pair,reference=p.reference)
    BasisPolicy(origin,rows,products,density,radial,mask,species,orders,pair,reference)
end
function model_with(m;radial=m.radial,pair_radial=m.pair_radial,radial_spec=m.radial_spec,
                    density=m.density,products=m.products,coupling=m.coupling,policy=m.policy)
    ResolvedModel(copy(m.atomic_numbers),copy(m.reference),radial,pair_radial,radial_spec,
        copy(m.pair_spec),m.lmax,density,products,coupling,m.source_sha256,policy)
end
basis_identity(m)=bytes2hex(SHA.sha256(string(m.policy.origin,m.policy.rows)))

"""Owning complete-invariant metadata. Correlation orders count density factors,
not distinct neighbors. `index_weight` is sum(n+l), not a physical polynomial
degree. Coupled terms retain their individual angular/radial dependencies.
"""
function descriptors(m)
    [begin
        columns=structural_columns(m,(i,))
        terms=[(;product=j,original_product=collect(m.policy.products)[j],coefficient=m.coupling[i,j],
            density=copy(m.products[j]),radial=[m.density[a][1] for a in m.products[j]],
            original_density=[m.policy.density[a] for a in m.products[j]],
            original_radial=[m.policy.radial[m.density[a][1]] for a in m.products[j]],
            harmonics=[m.density[a][2] for a in m.products[j]],
            angular=[(l=isqrt(m.density[a][2]-1),m=m.density[a][2]-(isqrt(m.density[a][2]-1)^2+isqrt(m.density[a][2]-1)+1)) for a in m.products[j]]) for j in columns]
        orders=sort!(unique([length(m.products[j]) for j in columns]))
        radial_weights=family_index_weights(m); known=!isnothing(radial_weights)
        weights=known ? [sum(radial_weights[m.density[a][1]] for a in m.products[j]) for j in columns] : Int[]
        (;id=(m.policy.origin,collect(m.policy.rows)[i]),index=i,original_index=collect(m.policy.rows)[i],
          enabled=collect(Bool,m.policy.mask)[i],orders,nominal_body_orders=orders .+ 1,
          index_weight=known ? maximum(weights;init=0) : nothing,
          normalization="real spherical L2; imported sparse coupling",terms,
          coefficient_indices=[i+(z-1)*size(m.coupling,1) for z in eachindex(m.atomic_numbers)],
          coefficient_ids=[(m.policy.origin,m.policy.rows[i],z) for z in m.atomic_numbers],
          species_orders=m.policy.orders,selected_species=collect(m.policy.species))
    end for i=1:size(m.coupling,1)]
end
function basis_summary(m,ps=nothing)
    ds=descriptors(m)
    counts=Dict{Int,Int}(); weights=Dict{Int,Int}()
    for d in ds
        for k in d.orders; counts[k]=get(counts,k,0)+1; end
        isnothing(d.index_weight) || (weights[d.index_weight]=get(weights,d.index_weight,0)+1)
    end
    (;identity=basis_identity(m),descriptors=ds,order_counts=counts,index_weight_counts=weights,
      provenance=(basis=m.policy.origin,invariants=collect(Int,m.policy.rows),products=collect(Int,m.policy.products),
                  density=collect(Int,m.policy.density),radial=collect(Int,m.policy.radial)),
      dimensions=(invariants=size(m.coupling,1),products=length(m.products),density=length(m.density),
                  radial=length(m.radial_spec),harmonics=(m.lmax+1)^2),
      coefficients=isnothing(ps) ? Float64[] : copy(vec(ps.WB)))
end
struct BasisSelection
    basis::String
    indices::Tuple
end
"""Select complete invariant rows by indices/IDs and optional order/index weight.
For mixed-order rows, `orders` requires every actual order to be in the criterion.
"""
function BasisSelection(m;indices=nothing,ids=nothing,orders=nothing,max_index_weight=nothing)
    ds=descriptors(m)
    isnothing(indices) || isnothing(ids) || throw(ArgumentError("choose indices or IDs"))
    rows=isnothing(indices) ? collect(1:length(ds)) : collect(indices)
    if !isnothing(ids)
        rows=Int[]
        for id in ids
            i=findfirst(d -> d.id==id,ds)
            isnothing(i) && throw(ArgumentError("descriptor ID belongs to another basis"))
            push!(rows,i)
        end
    end
    all(i -> i isa Integer && !(i isa Bool) && 1<=i<=length(ds),rows) && length(unique(rows))==length(rows) ||
        throw(ArgumentError("invalid/duplicate invariant indices"))
    if !isnothing(max_index_weight)
        any(d -> isnothing(d.index_weight),ds[rows]) && throw(ArgumentError("undefined radial index weight; select explicit IDs"))
    end
    filter!(i -> (isnothing(orders) || all(k -> k in orders,ds[i].orders)) &&
        (isnothing(max_index_weight) || ds[i].index_weight<=max_index_weight),rows)
    BasisSelection(basis_identity(m),Tuple(Int.(rows)))
end
function checked_selection(m,s)
    s.basis==basis_identity(m) || throw(ArgumentError("selection belongs to another basis"))
    all(i -> i isa Integer && !(i isa Bool) && 1<=i<=size(m.coupling,1),s.indices) && length(unique(s.indices))==length(s.indices) ||
        throw(ArgumentError("invalid selection"))
    collect(Int,s.indices)
end
"""Reversible effective coefficient mask; ps is never overwritten. Reapply an
all-row selection to remove the mask. Fitting columns and AD use the same mask.
"""
function mask_basis(m,s::BasisSelection)
    rows=checked_selection(m,s)
    model_with(m;policy=policy_with(m.policy;mask=Tuple(i in rows for i=1:size(m.coupling,1))))
end
function components(m;pair=m.policy.pair,reference=m.policy.reference)
    model_with(m;policy=policy_with(m.policy;pair,reference))
end
function project_species(m,species,orders::AbstractUnitRange{<:Integer})
    first(orders)>=0 && !isempty(orders) || throw(ArgumentError("orders must be a nonnegative inclusive range"))
    group=species isa Union{Symbol,ChemicalSpecies,Integer} ? (species,) : species
    channels=sort!(unique([element_channel(m,s) for s in group]))
    model_with(m;policy=policy_with(m.policy;species=Tuple(channels),orders=Int(first(orders)):Int(last(orders))))
end
project_species(m,species,orders)=throw(ArgumentError("orders must be a nonnegative inclusive unit range"))
order_allowed(m,k)=isnothing(m.policy.orders) || k in m.policy.orders
pair_allowed(m,z)=m.policy.pair && order_allowed(m,z in m.policy.species ? 1 : 0)
reference_value(m,z)=m.policy.reference && order_allowed(m,0) ? m.reference[z] : zero(m.reference[z])
effective_WB(m,ps)=ps.WB .* collect(Bool,m.policy.mask)

"""Structural selection from an existing basis. Returns model/parameters and
explicit current-to-parent dependency maps. No evaluator/workspace is constructed.
Pair and references are preserved unless component flags were explicitly changed.
"""
function prune_basis(m,ps,s::BasisSelection)
    rows=checked_selection(m,s)
    cols=structural_columns(m,rows)
    den=sort!(unique(vcat(Int[],m.products[cols]...)))
    rad=sort!(unique([m.density[a][1] for a in den]))
    dmap=Dict(a=>i for (i,a) in enumerate(den)); rmap=Dict(a=>i for (i,a) in enumerate(rad))
    channels=map(ch -> subset_radial(ch,rad),m.radial)
    policy=policy_with(m.policy;rows=Tuple(collect(m.policy.rows)[rows]),products=Tuple(collect(m.policy.products)[cols]),
        density=Tuple(collect(m.policy.density)[den]),radial=Tuple(collect(m.policy.radial)[rad]),mask=Tuple(collect(Bool,m.policy.mask)[rows]))
    out=model_with(m;radial=channels,radial_spec=m.radial_spec[rad],
        density=[(rmap[m.density[a][1]],m.density[a][2]) for a in den],
        products=[[dmap[a] for a in m.products[j]] for j in cols],coupling=m.coupling[rows,cols],policy)
    params=merge(deepcopy(ps),(WB=copy(ps.WB[rows,:]),))
    hasproperty(ps,:radial) && (params=merge(params,(radial=merge(ps.radial,(scale=ps.radial.scale[rad,:,:],)),)))
    maps=(invariants=rows,products=cols,density=den,radial=rad,harmonics=sort!(unique(last.(out.density))))
    (;model=out,parameters=params,maps,coefficient_map=coefficient_map(m,rows))
end
restrict_basis(m,ps,s::BasisSelection)=prune_basis(m,ps,s)
function restrict_basis(path::AbstractString;kwargs...)
    # Inspect the serialized basis and compact its numerical tables BEFORE
    # constructing ResolvedModel. No full radial objects, evaluator or executor.
    d=TOML.parsefile(path)
    get(d,"schema",nothing)=="PotentialFrameWork.resolved_ACE" && get(d,"schema_version",nothing)==1 ||
        throw(ArgumentError("restricted construction requires a resolved v1 basis specification"))
    d["radial"]["kind"]=="piecewise_cubic_agnesi" || throw(ArgumentError("use restrict_basis(model,ps,selection) for custom bases"))
    raw=d["coupling"]; shape=Int.(raw["shape"])
    C=sparse(Int.(raw["rows"]),Int.(raw["columns"]),Float64.(raw["values"]),shape...)
    rspec=[(Int(s["n"]),Int(s["l"])) for s in d["radial"]["spec"]]
    density=Tuple{Int,Int}.(d["density_indices"]); products=[Int.(v) for v in d["products"]]
    policy=default_policy(d["source_sha256"],rspec,density,products,C)
    spec=(;coupling=C,radial_spec=rspec,density,products,policy,
        atomic_numbers=Int.(d["species_atomic_numbers"]),pair_spec=[(Int(v["n"]),Int(v["l"])) for v in d["pair_radial"]["spec"]],
        radial=[SpecificationRadial(sum.(rspec))])
    haskey(d,"interface_policy") && (spec=merge(spec,(policy=read_policy(d["interface_policy"],spec),)))
    rows=checked_selection(spec,BasisSelection(spec;kwargs...)); cols=structural_columns(spec,rows)
    den=sort!(unique(vcat(Int[],products[cols]...))); rad=sort!(unique(first.(density[den])))
    dmap=Dict(a=>i for (i,a) in enumerate(den)); rmap=Dict(a=>i for (i,a) in enumerate(rad))
    for ch in d["radial"]["channels"]
        coeff=reshape(Float64.(ch["coefficients"]),Tuple(Int.(ch["coefficient_shape"])))[:,rad,:]
        ch["coefficient_shape"]=collect(size(coeff)); ch["coefficients"]=vec(coeff)
    end
    d["radial"]["spec"]=d["radial"]["spec"][rad]
    d["density_indices"]=[[rmap[density[a][1]],density[a][2]] for a in den]
    d["products"]=[[dmap[a] for a in products[j]] for j in cols]
    reduced=C[rows,cols]; ri,ci,vals=findnz(reduced)
    d["coupling"]=Dict("shape"=>collect(size(reduced)),"rows"=>ri,"columns"=>ci,"values"=>vals)
    params=d["parameters"]; oldshape=Tuple(Int.(params["WB_shape"]))
    count=prod(oldshape); flat=Float64.(params["flat"])
    WB=reshape(flat[1:count],oldshape)[rows,:]
    params["WB_shape"]=collect(size(WB)); params["flat"]=vcat(vec(WB),flat[count+1:end])
    p=spec.policy
    d["interface_policy"]=checkpoint_policy(policy_with(p;rows=Tuple(collect(p.rows)[rows]),
        products=Tuple(collect(p.products)[cols]),density=Tuple(collect(p.density)[den]),
        radial=Tuple(collect(p.radial)[rad]),mask=Tuple(collect(Bool,p.mask)[rows])))
    model,parameters,state=read_resolved(d)
    maps=(invariants=rows,products=cols,density=den,radial=rad,harmonics=sort!(unique(last.(model.density))))
    (;model,parameters,state,maps,coefficient_map=coefficient_map(spec,rows))
end
struct SpecificationRadial
    weights::Vector{Int}
end
radial_index_weights(ch::SpecificationRadial)=ch.weights
function coefficient_map(m,rows;fit_pair=true)
    nb=size(m.coupling,1); nz=length(m.atomic_numbers); nq=length(m.pair_spec)
    indices=vcat([i+(z-1)*nb for z=1:nz for i in rows],fit_pair ? collect(nb*nz+1:(nb+nq)*nz) : Int[])
    sparse(indices,collect(eachindex(indices)),ones(length(indices)),(nb+(fit_pair ? nq : 0))*nz,length(indices))
end
mapped_prior(R,T)=R*T
"""Differentiable numerical gather into a pruned basis's parent maps. Fixed
structural maps stay discrete; Zygote scatters readout/radial-scale adjoints back
to parent coordinates. Unchanged subtrees are shared; this function does not
mutate them. `prune_basis` itself returns owning parameters instead.
"""
function restrict_parameters(ps,maps)
    result=merge(ps,(WB=ps.WB[maps.invariants,:],))
    if hasproperty(ps,:radial)
        result=merge(result,(radial=merge(ps.radial,(scale=ps.radial.scale[maps.radial,:,:],)),))
    end
    result
end
function center_indices(centers,n)
    isnothing(centers) && return 1:n
    c=collect(centers)
    all(i -> i isa Integer && !(i isa Bool) && 1<=i<=n,c) && length(unique(c))==length(c) ||
        throw(ArgumentError("centers must be unique original atom indices in 1:$n"))
    Int.(c)
end

function plot_basis(args...;kwargs...)
    throw(ArgumentError("load UnicodePlots for ACE plots"))
end

function structural_columns(m,rows)
    [j for j=1:size(m.coupling,2) if any(k -> rowvals(m.coupling)[k] in rows,nzrange(m.coupling,j))]
end
