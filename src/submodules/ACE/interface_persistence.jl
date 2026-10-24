function checkpoint_policy(p)
    Dict("version"=>1,"origin"=>p.origin,"rows"=>collect(Int,p.rows),
        "products"=>collect(Int,p.products),"density"=>collect(Int,p.density),
        "radial"=>collect(Int,p.radial),"mask"=>collect(Bool,p.mask),"species"=>collect(Int,p.species),
        "orders"=>isnothing(p.orders) ? Int[] : [first(p.orders),last(p.orders)],
        "pair"=>p.pair,"reference"=>p.reference)
end
function read_policy(d,m)
    get(d,"version",nothing)==1 || throw(ArgumentError("unsupported basis policy version"))
    for (key,n) in (("rows",size(m.coupling,1)),("products",length(m.products)),
                    ("density",length(m.density)),("radial",length(m.radial_spec)))
        v=d[key]
        length(v)==n && all(x -> x isa Integer && x>0,v) && length(unique(v))==n ||
            throw(ArgumentError("invalid $key provenance map"))
    end
    length(d["mask"])==size(m.coupling,1) && all(x -> x isa Bool,d["mask"]) || throw(ArgumentError("invalid invariant mask"))
    sp=Int.(d["species"])
    all(i -> 1<=i<=length(m.atomic_numbers),sp) && length(unique(sp))==length(sp) || throw(ArgumentError("invalid projected species"))
    o=d["orders"]
    isempty(o) || (length(o)==2 && 0<=o[1]<=o[2]) || throw(ArgumentError("invalid projected orders"))
    length(d["origin"])==64 || throw(ArgumentError("invalid basis namespace"))
    BasisPolicy(d["origin"],Tuple(Int.(d["rows"])),Tuple(Int.(d["products"])),Tuple(Int.(d["density"])),
        Tuple(Int.(d["radial"])),Tuple(Bool.(d["mask"])),Tuple(sp),isempty(o) ? nothing : (Int(o[1]):Int(o[2])),
        Bool(d["pair"]),Bool(d["reference"]))
end
function checkpoint_metadata(m::RadialMetadata)
    Dict("spec"=>collect.(m.spec),"cutoff"=>m.cutoff,
        "index_weights"=>isnothing(m.index_weights) ? Int[] : m.index_weights,
        "weights_known"=>!isnothing(m.index_weights),"simd"=>m.simd,"ad"=>m.ad,"origin"=>string(m.origin),
        "envelope"=>string(m.envelope),"channel_dependence"=>string(m.channel_dependence))
end
function read_metadata(d)
    RadialMetadata(Tuple{Int,Int}.(d["spec"]);cutoff=d["cutoff"],
        index_weights=d["weights_known"] ? d["index_weights"] : nothing,
        simd=d["simd"],ad=d["ad"],origin=Symbol(d["origin"]),
        envelope=Symbol(d["envelope"]),channel_dependence=Symbol(d["channel_dependence"]))
end
function checkpoint_channel(ch::CallableRadial,ps)
    isnothing(ch.persistence) && throw(ArgumentError("unregistered callable: provide a versioned registration/source bundle or explicitly tabulate it"))
    record=ch.persistence; route=get(record,:route,"")
    d=Dict{String,Any}("kind"=>"callable","metadata"=>checkpoint_metadata(ch.metadata),
        "config"=>checkpoint_encode(ch.config),"route"=>route,"interface_version"=>1)
    if route=="registered"
        haskey(RADIAL_REGISTRY,record.id) || throw(ArgumentError("unregistered radial ID $(record.id)"))
        d["id"]=record.id
    elseif route=="source"
        record.interface_version==1 || throw(ArgumentError("unsupported source interface"))
        bytes2hex(SHA.sha256(record.source))==record.source_sha256 || throw(ArgumentError("custom source hash mismatch"))
        d["source"]=record.source; d["entrypoint"]=record.entrypoint
        d["source_sha256"]=record.source_sha256; d["dependencies"]=record.dependencies
    else
        throw(ArgumentError("unknown callable persistence route"))
    end
    samples=map((0.23,0.59,0.83)) do fraction
        r=fraction*radial_cutoff(ch)
        v,dr=Base.invokelatest(radial_ed,ch,r,ps)
        all(isfinite,v) && all(isfinite,dr) || throw(ArgumentError("nonfinite custom radial reference"))
        Dict("r"=>r,"values"=>v,"derivatives"=>dr)
    end
    d["reference"]=collect(samples)
    d
end
function checkpoint_channel(ch::SubsetRadial,ps)
    Dict("kind"=>"subset","parent"=>checkpoint_channel(ch.parent,ps),"indices"=>ch.indices)
end
function checkpoint_channel(ch::TabulatedRadial,ps)
    Dict("kind"=>"raw_distance_cubic_v1","metadata"=>checkpoint_metadata(ch.metadata),
        "nodes"=>ch.nodes,"coefficients"=>checkpoint_encode(ch.coefficients),"report"=>checkpoint_encode(ch.report))
end
function checkpoint_radial(channels,spec;parameters=NamedTuple())
    Dict("kind"=>"interface_channels_v1","spec"=>[Dict("n"=>n,"l"=>l) for (n,l) in spec],
        "shape"=>collect(size(channels)),"channels"=>[checkpoint_channel(ch,parameters) for ch in vec(channels)])
end
function read_channel(d;allow_code=false,parameters=NamedTuple(),source_cache=Dict())
    kind=d["kind"]
    if kind=="subset"
        parent=read_channel(d["parent"];allow_code,parameters,source_cache); indices=Int.(d["indices"])
        all(i -> 1<=i<=length(channel_spec(parent)),indices) && length(unique(indices))==length(indices) || throw(ArgumentError("invalid radial subset"))
        return subset_radial(parent,indices)
    elseif kind=="raw_distance_cubic_v1"
        metadata=read_metadata(d["metadata"]); xs=Float64.(d["nodes"])
        length(xs)>=2 && first(xs)>=0 && last(xs)==metadata.cutoff && all(>(0),diff(xs)) || throw(ArgumentError("invalid tabulation nodes"))
        c=checkpoint_decode(d["coefficients"])
        size(c)==(4,length(metadata.spec),length(xs)-1) || throw(DimensionMismatch("tabulation coefficients"))
        return TabulatedRadial(xs,c,metadata,checkpoint_decode(d["report"]))
    elseif kind!="callable"
        throw(ArgumentError("unsupported radial channel kind $kind"))
    end
    d["interface_version"]==1 || throw(ArgumentError("unsupported custom radial interface"))
    # Refusal precedes parsing/evaluating ANY saved source.
    d["route"]=="source" && !allow_code && throw(ArgumentError("custom Julia source loading is disabled; pass allow_code=true for trusted code"))
    metadata=read_metadata(d["metadata"]); config=checkpoint_decode(d["config"])
    ch=if d["route"]=="registered"
        factory=get(RADIAL_REGISTRY,d["id"],nothing)
        isnothing(factory) && throw(ArgumentError("radial implementation $(d["id"]) is not registered; load its extension first"))
        CallableRadial(factory(config,metadata),metadata;config,persistence=(route="registered",id=d["id"]))
    elseif d["route"]=="source"
        bytes2hex(SHA.sha256(d["source"]))==d["source_sha256"] || throw(ArgumentError("custom source hash mismatch"))
        source_radial(d["source"],d["entrypoint"],metadata;config,dependencies=Dict{String,String}(d["dependencies"]),allow_code,namespace_cache=source_cache)
    else
        throw(ArgumentError("unknown custom radial route"))
    end
    for ref in d["reference"]
        v,dr=Base.invokelatest(radial_ed,ch,ref["r"],parameters)
        length(v)==length(ref["values"]) && length(dr)==length(ref["derivatives"]) || throw(DimensionMismatch("custom radial reference dimensions"))
        isapprox(v,ref["values"];atol=5e-10,rtol=5e-9) && isapprox(dr,ref["derivatives"];atol=5e-10,rtol=5e-9) ||
            throw(ArgumentError("custom radial numerical reference mismatch"))
    end
    caps=Base.invokelatest(radial_capabilities,ch,parameters)
    caps.cpu_scalar && (!metadata.simd || caps.cpu_simd) || throw(ArgumentError("custom radial capability validation failed: $(caps.reasons)"))
    metadata.ad && caps.training_ad_checked && !caps.training_ad_local &&
        throw(ArgumentError("custom radial declared AD capability failed: $(caps.reasons)"))
    ch
end
function read_custom_radial(d,nz;allow_code=false,parameters=NamedTuple())
    Tuple(d["shape"])==(nz,nz) && length(d["channels"])==nz^2 || throw(DimensionMismatch("custom species channels"))
    spec=[(Int(v["n"]),Int(v["l"])) for v in d["spec"]]
    source_cache=Dict()
    channels=[read_channel(ch;allow_code,parameters,source_cache) for ch in d["channels"]]
    all(ch -> channel_spec(ch)==spec,channels) || throw(ArgumentError("custom radial spec mismatch"))
    # Homogeneous channels retain concrete hot-loop types; heterogeneous callable
    # families must be adapted to one concretely typed implementation.
    length(unique(typeof.(channels)))==1 || throw(ArgumentError("custom channel types must be homogeneous"))
    reshape(typeof(first(channels))[channels...],nz,nz),spec
end
channel_spec(ch)=ch.metadata.spec
channel_spec(ch::SubsetRadial)=channel_spec(ch.parent)[ch.indices]
function validate_checkpoint_model(m,ps)
    nz=length(m.atomic_numbers)
    length(unique(m.atomic_numbers))==nz && all(>(0),m.atomic_numbers) || throw(ArgumentError("invalid species channels"))
    length(m.reference)==nz && all(isfinite,m.reference) && all(isfinite,nonzeros(m.coupling)) || throw(ArgumentError("nonfinite/invalid model values"))
    size(m.radial)==size(m.pair_radial)==(nz,nz) || throw(DimensionMismatch("radial species channels"))
    size(m.coupling,2)==length(m.products) || throw(DimensionMismatch("coupling products"))
    all(indices -> !isempty(indices) && all(i -> 1<=i<=length(m.density),indices),m.products) || throw(ArgumentError("invalid product indices"))
    all(v -> 1<=v[1]<=length(m.radial_spec) && 1<=v[2]<=(m.lmax+1)^2,m.density) || throw(ArgumentError("invalid density indices"))
    read_policy(checkpoint_policy(m.policy),m)
    size(ps.WB)==(size(m.coupling,1),nz) && size(ps.Wpair)==(length(m.pair_spec),nz) || throw(DimensionMismatch("checkpoint readout dimensions"))
    for (key,n) in ((:radial,length(m.radial_spec)),(:pair,length(m.pair_spec)))
        hasproperty(ps,key) || continue
        size(getproperty(ps,key).scale)==(n,nz,nz) || throw(DimensionMismatch("checkpoint $key scale dimensions"))
    end
    nothing
end
