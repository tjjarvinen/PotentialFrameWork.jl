struct AgnesiTransform{T}
    p::Int
    q::Int
    a::T
    rin::T
    r0::T
    yin::T
    ycut::T
    rcut::T
end
struct RadialChannel{T}
    transform::AgnesiTransform{T}
    envelope::Symbol
    envelope_parameters::Vector{T}
    nodes::Vector{T}
    coefficients::Array{T,3} # polynomial power, radial channel, interval
end
struct ResolvedModel{T}
    atomic_numbers::Vector{Int}
    reference::Vector{T}
    radial::Matrix{RadialChannel{T}}
    pair_radial::Matrix{RadialChannel{T}}
    radial_spec::Vector{Tuple{Int,Int}}
    pair_spec::Vector{Tuple{Int,Int}}
    lmax::Int
    density::Vector{Tuple{Int,Int}}
    products::Vector{Vector{Int}}
    coupling::SparseMatrixCSC{T,Int}
    source_sha256::String
end

function read_radial(D,nz)
    D["kind"]=="piecewise_cubic_agnesi" || throw(ArgumentError("unsupported resolved radial kind"))
    result=Matrix{RadialChannel{Float64}}(undef,nz,nz)
    seen=Set{Tuple{Int,Int}}()
    spec=[(Int(b["n"]),Int(b["l"])) for b in D["spec"]]
    for ch in D["channels"]
        ij=(Int(ch["centre"]),Int(ch["neighbour"]))
        all(i -> 1<=i<=nz,ij) && !(ij in seen) || throw(ArgumentError("invalid/duplicate radial species channel"))
        push!(seen,ij)
        t=ch["transform"]
        t["normalization_rin"]==t["rin"] || throw(ArgumentError("unsupported differing transform origins"))
        transform=AgnesiTransform(Int(t["p"]),Int(t["q"]),t["a"],t["rin"],t["r0"],t["yin"],t["ycut"],t["rcut"])
        transform.p>0 && transform.q>=transform.p && transform.a>0 && transform.rcut>transform.r0>transform.rin ||
            throw(ArgumentError("invalid Agnesi transform"))
        kind=Symbol(ch["envelope_kind"]); e=ch["envelope"]
        ep=if kind===:PolyEnvelope2sX
            [e["x1"],e["x2"],e["p1"],e["p2"],e["s"]]
        elseif kind===:ACE1_PolyEnvelope1sR
            [e["rcut"],e["r0"],e["p"]]
        else
            throw(ArgumentError("unsupported resolved envelope $kind"))
        end
        nodes=Float64.(ch["nodes"]); shape=Tuple(Int.(ch["coefficient_shape"]))
        length(nodes)>=2 && first(nodes)==-1 && last(nodes)==1 && issorted(nodes) && all(>(0),diff(nodes)) ||
            throw(ArgumentError("invalid spline nodes"))
        shape==(4,length(spec),length(nodes)-1) || throw(DimensionMismatch("invalid spline coefficients"))
        coeff=reshape(Float64.(ch["coefficients"]),shape)
        all(isfinite,coeff) || throw(ArgumentError("nonfinite spline coefficients"))
        result[ij...]=RadialChannel(transform,kind,Float64.(ep),nodes,coeff)
    end
    length(seen)==nz^2 || throw(ArgumentError("missing radial species channels"))
    result,spec
end

"""Load only the versioned resolved schema, without importing ACEpotentials.
Return `(model, parameters, state)`; parameters own the original readout matrices.
The recipe JSON is exported separately in its pinned producer environment.
"""
function load_resolved(path)
    D=TOML.parsefile(path)
    get(D,"schema",nothing)=="PotentialFrameWork.resolved_ACE" && get(D,"schema_version",nothing)==1 ||
        throw(ArgumentError("unsupported resolved ACE schema; use the pinned export_resolved_model.jl producer"))
    z=Int.(D["species_atomic_numbers"]); nz=length(z)
    length(unique(z))==nz && all(>(0),z) || throw(ArgumentError("invalid species channels"))
    angular=D["angular"]
    angular["kind"]=="real_spherical" && angular["normalization"]=="L2" && angular["ordering"]=="l*l+l+m+1" ||
        throw(ArgumentError("unsupported angular convention"))
    lmax=Int(angular["lmax"]); lmax>=0 || throw(ArgumentError("invalid lmax"))
    radial,rspec=read_radial(D["radial"],nz); pair,pspec=read_radial(D["pair_radial"],nz)
    density=[Tuple(Int.(v)) for v in D["density_indices"]]
    all(v -> length(v)==2 && 1<=v[1]<=length(rspec) && 1<=v[2]<=(lmax+1)^2,density) ||
        throw(ArgumentError("invalid density indices"))
    products=[Int.(v) for v in D["products"]]
    all(v -> !isempty(v) && all(i -> 1<=i<=length(density),v),products) ||
        throw(ArgumentError("invalid product indices"))
    c=D["coupling"]; m,n=Int.(c["shape"])
    n==length(products) || throw(DimensionMismatch("coupling/product dimensions"))
    C=sparse(Int.(c["rows"]),Int.(c["columns"]),Float64.(c["values"]),m,n)
    params=D["parameters"]; flat=Float64.(params["flat"])
    params["ordering"]=="vec(WB), vec(Wpair); Julia column-major" || throw(ArgumentError("unsupported parameter ordering"))
    params["WB_shape"]==[m,nz] && params["Wpair_shape"]==[length(pspec),nz] ||
        throw(DimensionMismatch("readout dimensions"))
    length(flat)==(m+length(pspec))*nz || throw(DimensionMismatch("parameter length"))
    length(D["reference_eV"])==nz || throw(DimensionMismatch("reference energy channels"))
    all(isfinite,flat) && all(isfinite,D["reference_eV"]) && all(isfinite,nonzeros(C)) || throw(ArgumentError("nonfinite model values"))
    WB=reshape(copy(flat[1:m*nz]),m,nz)
    Wpair=reshape(copy(flat[m*nz+1:end]),length(pspec),nz)
    model=ResolvedModel(z,Float64.(D["reference_eV"]),radial,pair,rspec,pspec,lmax,density,products,C,D["source_sha256"])
    model,(;WB,Wpair),NamedTuple()
end

"""Convenience import for the supplied recipe with its separately produced
resolved export. The source hash must match; otherwise rerun the pinned exporter.
No ACEpotentials constructor or evaluator is loaded in the package process.
"""
function load_recipe(path;resolved_path=joinpath(dirname(path),"resolved_model_v1.toml"))
    isfile(resolved_path) || throw(ArgumentError("missing resolved export; run data/ace_ref/export_resolved_model.jl in its pinned environment"))
    model,ps,st=load_resolved(resolved_path)
    bytes2hex(SHA.sha256(read(path)))==model.source_sha256 || throw(ArgumentError("recipe hash differs from the resolved export; regenerate it in the pinned environment"))
    model,ps,st
end

function transformed(t::AgnesiTransform,r)
    s=(r-t.rin)/(t.r0-t.rin)
    y=r<=t.rin ? one(r) : inv(1+t.a*s^t.q/(1+s^(t.q-t.p)))
    clamp(-1+2*(y-t.yin)/(t.ycut-t.yin),-one(y),one(y))
end
function envelope(ch,r,x)
    e=ch.envelope_parameters
    if ch.envelope===:PolyEnvelope2sX
        e[1]<x<e[2] || return zero(x)
        e[5]*(x-e[1])^Int(e[3])*(e[2]-x)^Int(e[4])
    else
        r>e[1] && return zero(r)
        r>0 || throw(DomainError(r,"singular imported pair envelope at coincident atoms"))
        s=r/e[2]; sc=e[1]/e[2]; p=Int(e[3])
        s^(-p)-sc^(-p)+p*sc^(-p-1)*(s-sc)
    end
end

"""Independent scalar evaluation of the resolved spline and envelope.
Preserves promoted arithmetic (including dual numbers); numeric Å/eV conventions.
"""
function radial_values(ch::RadialChannel,r::Real)
    x=transformed(ch.transform,r)
    k=clamp(searchsortedlast(ch.nodes,x),1,length(ch.nodes)-1)
    t=(x-ch.nodes[k])/(ch.nodes[k+1]-ch.nodes[k])
    e=envelope(ch,r,x); c=ch.coefficients
    [(((c[4,n,k]*t+c[3,n,k])*t+c[2,n,k])*t+c[1,n,k])*e for n=1:size(c,2)]
end
