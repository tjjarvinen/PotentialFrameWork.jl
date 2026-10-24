"""Explicit channel order and capabilities. Callables own their envelope once.
`spec` lists (n,l) labels; optional index weights are algebraic prior metadata.
Cutoff is fixed in Å; `origin` is a declared policy, not an inferred dimension.
"""
struct RadialMetadata
    spec::Vector{Tuple{Int,Int}}
    cutoff::Float64
    index_weights::Union{Nothing,Vector{Int}}
    simd::Bool
    ad::Bool
    origin::Symbol
    envelope::Symbol
    channel_dependence::Symbol
end
function RadialMetadata(spec;cutoff,index_weights=nothing,simd=false,ad=false,origin=:error,
                        envelope=:included,channel_dependence=:central_neighbor)
    cutoff>0 && isfinite(cutoff) || throw(ArgumentError("positive finite radial cutoff required"))
    origin in (:error,:regular) || throw(ArgumentError("origin policy must be :error or :regular"))
    envelope===:included || throw(ArgumentError("callables must include their envelope exactly once"))
    channel_dependence in (:central_neighbor,:independent) || throw(ArgumentError("unsupported channel dependence"))
    sp=Tuple{Int,Int}.(spec)
    iw=isnothing(index_weights) ? nothing : Int.(index_weights)
    isnothing(iw) || length(iw)==length(sp) || throw(DimensionMismatch("radial index weights"))
    RadialMetadata(sp,Float64(cutoff),iw,simd,ad,origin,envelope,channel_dependence)
end
struct CallableRadial{F,C,S}
    evaluator::F
    config::C
    metadata::RadialMetadata
    persistence::S
end
# Code closures commonly capture a module binding (e.g. a local `A=ACE`). Julia
# cannot deepcopy modules. Share immutable code references while still owning
# every mutable config/metadata/parameter value saved by a pullback. Stateful
# evaluators that can be copied are copied; noncopyable mutable captures fail
# explicitly instead of borrowing mutable forward inputs.
function immutable_radial_code(x)
    x isa Union{Module,Type,Symbol,String} && return true
    isbitstype(typeof(x)) && return true
    ismutabletype(typeof(x)) && return false
    all(i -> immutable_radial_code(getfield(x,i)),1:fieldcount(typeof(x)))
end
function Base.deepcopy_internal(ch::CallableRadial,dict::IdDict)
    haskey(dict,ch) && return dict[ch]
    evaluator=if immutable_radial_code(ch.evaluator)
        ch.evaluator
    else
        try
            Base.deepcopy_internal(ch.evaluator,dict)
        catch e
            throw(ArgumentError("cannot snapshot radial callable; move mutable captures into explicit config/ps: "*sprint(showerror,e)))
        end
    end
    result=CallableRadial(evaluator,Base.deepcopy_internal(ch.config,dict),
        Base.deepcopy_internal(ch.metadata,dict),Base.deepcopy_internal(ch.persistence,dict))
    dict[ch]=result
    result
end
"""Wrap an in-place callable `(values, derivatives_or_nothing, r, ps, config)`.
Outputs follow metadata.spec. Generic code may accept independent SIMD distances;
the CPU adapter currently evaluates callable channels lane-wise (no SIMD AD
promise). Use `derivative=:ad` to adapt a values-only callable with ForwardDiff.
"""
function CallableRadial(f,metadata::RadialMetadata;config=NamedTuple(),persistence=nothing,derivative=:analytic)
    derivative in (:analytic,:ad) || throw(ArgumentError("unknown derivative route"))
    evaluator=derivative===:ad ? ADRadial(f) : f
    CallableRadial(evaluator,deepcopy(config),metadata,deepcopy(persistence))
end
struct ADRadial{F}
    values!::F
end
radial_ad!(args...)=throw(ArgumentError("the scalar radial AD fallback requires ForwardDiff and a scalar distance; use scalar lanes for SIMD inputs"))
(f::ADRadial)(v,d,r,ps,c)=isnothing(d) ? f.values!(v,r,ps,c) : radial_ad!(v,d,f.values!,r,ps,c)
function radial_evaluate!(v,d,ch::CallableRadial,r,ps=NamedTuple())
    length(v)==length(ch.metadata.spec) && (isnothing(d) || length(d)==length(v)) || throw(DimensionMismatch("radial output shape"))
    ch.evaluator(v,d,promoted_distance(r,radial_parameter_seed(ps)),ps,ch.config)
    length(v)==length(ch.metadata.spec) && (isnothing(d) || length(d)==length(v)) || throw(DimensionMismatch("callable changed radial output dimensions"))
    nothing
end
promoted_distance(r,seed)=r+seed
function promoted_distance(r::Vec{W},seed) where W
    seed isa Union{Integer,Float32,Float64} || throw(ArgumentError("custom SIMD AD is unsupported; use scalar lanes"))
    Vec(ntuple(i -> r[i]+seed,W))
end
radial_cutoff(ch::RadialChannel)=ch.transform.rcut
radial_cutoff(ch::CallableRadial)=ch.metadata.cutoff
radial_index_weights(ch::RadialChannel)=nothing # supplied by model.radial_spec below
radial_index_weights(ch::CallableRadial)=ch.metadata.index_weights
function radial_ed(ch::CallableRadial,r::Real,ps=NamedTuple())
    T=typeof(r+radial_parameter_seed(ps)); n=length(ch.metadata.spec)
    v=Vector{T}(undef,n); d=similar(v)
    radial_evaluate!(v,d,ch,r,ps)
    v,d
end
function radial_values(ch::CallableRadial,r::Real,ps=NamedTuple())
    T=typeof(r+radial_parameter_seed(ps)); v=Vector{T}(undef,length(ch.metadata.spec))
    radial_evaluate!(v,nothing,ch,r,ps); v
end
radial_ed(ch::RadialChannel,r,ps)=radial_ed(ch,r)
radial_values(ch::RadialChannel,r,ps)=radial_values(ch,r)
radial_parameters(ps,key)=hasproperty(ps,key) ? getproperty(ps,key) : NamedTuple()
function family_index_weights(model;pair=false)
    spec=pair ? model.pair_spec : model.radial_spec
    channels=pair ? model.pair_radial : model.radial
    weights=nothing
    for ch in channels
        current=ch isa RadialChannel ? sum.(spec) : radial_index_weights(ch)
        isnothing(current) && return nothing
        length(current)==length(spec) || throw(DimensionMismatch("radial index weights"))
        if isnothing(weights)
            weights=copy(current)
        elseif weights!=current
            return nothing # no common declared index weight across channels
        end
    end
    weights
end

struct SubsetRadial{R}
    parent::R
    indices::Vector{Int}
end
subset_radial(ch::RadialChannel,indices)=RadialChannel(ch.transform,ch.envelope,copy(ch.envelope_parameters),copy(ch.nodes),ch.coefficients[:,indices,:])
subset_radial(ch,indices)=SubsetRadial(ch,copy(indices))
radial_cutoff(ch::SubsetRadial)=radial_cutoff(ch.parent)
radial_index_weights(ch::SubsetRadial)=isnothing(radial_index_weights(ch.parent)) ? nothing : radial_index_weights(ch.parent)[ch.indices]
function radial_ed(ch::SubsetRadial,r,ps=NamedTuple())
    v,d=radial_ed(ch.parent,r,ps); v[ch.indices],d[ch.indices]
end
radial_values(ch::SubsetRadial,r,ps=NamedTuple())=radial_values(ch.parent,r,ps)[ch.indices]

function replace_radial(m,ps,channels;parameters=NamedTuple(),pair=false)
    size(channels)==size(m.radial) || throw(DimensionMismatch("radial species channels"))
    spec=pair ? m.pair_spec : m.radial_spec
    all(ch -> channel_spec(ch)==spec,channels) || throw(ArgumentError("custom radial ordering must match existing (n,l) specification"))
    key=pair ? :custom_pair : :custom_radial
    identity=bytes2hex(SHA.sha256(string(m.policy.origin,pair,[radial_definition_id(ch) for ch in channels])))
    policy=policy_with(m.policy;origin=identity)
    model=pair ? model_with(m;pair_radial=copy(channels),policy) : model_with(m;radial=copy(channels),policy)
    (;model,parameters=merge(deepcopy(ps),NamedTuple{(key,)}((deepcopy(parameters),))))
end

"""Separate declared capabilities and diagnostic checks; GPU is always unverified.
This probes values/first derivatives; AD/force-loss requires generic evaluator
code and the independent derivative-product checks documented in the guide.
"""
function radial_capabilities(ch::CallableRadial,ps=NamedTuple();r=0.37radial_cutoff(ch))
    reasons=String[]; scalar=false; simd=false
    try
        v,d=radial_ed(ch,r,ps); scalar=all(isfinite,v) && all(isfinite,d)
        scalar || push!(reasons,"nonfinite scalar values/derivatives")
    catch e
        push!(reasons,sprint(showerror,e))
    end
    if ch.metadata.simd
        try
            v=fill(Vec((r,r+0.01)),length(ch.metadata.spec)); d=copy(v)
            radial_evaluate!(v,d,ch,Vec((r,r+0.01)),ps)
            simd=all(x -> all(isfinite,Tuple(x)),v) && all(x -> all(isfinite,Tuple(x)),d)
            for lane=1:2
                sv,sd=radial_ed(ch,r+(lane-1)*0.01,ps)
                simd &= isapprox([x[lane] for x in v],sv;atol=1e-10,rtol=1e-8) &&
                        isapprox([x[lane] for x in d],sd;atol=1e-10,rtol=1e-8)
            end
            simd || push!(reasons,"SIMD values/derivatives disagree with scalar samples")
        catch e
            push!(reasons,"SIMD: "*sprint(showerror,e))
        end
    else
        push!(reasons,"SIMD not declared; CPU uses scalar lanes")
    end
    ad=ch.metadata.ad ? radial_ad_probe(ch,ps,r) : (checked=false,passed=false,reason="training AD not declared")
    isempty(ad.reason) || push!(reasons,ad.reason)
    (;cpu_scalar=scalar,cpu_simd=simd,training_ad_declared=ch.metadata.ad,
      training_ad_checked=ad.checked,training_ad_local=ad.passed,gpu_verified=false,reasons)
end
function radial_samples(m;family=:radial,centre=1,neighbour=1,channels=nothing,points=100,parameters=NamedTuple())
    family in (:radial,:pair_radial) || throw(ArgumentError("unknown radial family"))
    ch=getproperty(m,family)[centre,neighbour]
    channels=isnothing(channels) ? (1:min(6,length(family===:radial ? m.radial_spec : m.pair_spec))) : channels
    points>=2 || throw(ArgumentError("need at least two points"))
    x=collect(range(0.01radial_cutoff(ch),radial_cutoff(ch);length=points))
    (;distances=x,values=permutedims(hcat([radial_values(ch,r,parameters)[channels] for r in x]...)),channels=collect(channels))
end

# A portable, smooth reference callable used in examples and registration tests.
struct SmoothExponential end
function radial_power(x,n)
    value=one(x)
    for _=1:n; value*=x; end
    value
end
function (::SmoothExponential)(v,d,r,ps,c)
    rc=c.cutoff; x=r/rc; q=max(one(r)-x,zero(r))
    rate=ps.rate[1]; amplitude=ps.amplitude[1]; ex=exp(-rate*r)
    q2=q*q; q3=q2*q; q4=q2*q2
    for i in eachindex(v)
        n=c.powers[i]; xn=radial_power(x,n)
        v[i]=amplitude*xn*q4*ex
        if !isnothing(d)
            dx=n==0 ? zero(r) : n*radial_power(x,n-1)/rc
            d[i]=amplitude*ex*(dx*q4-4xn*q3/rc-rate*xn*q4)
        end
    end
    nothing
end

struct TabulatedRadial{T}
    nodes::Vector{T}
    coefficients::Array{T,3}
    metadata::RadialMetadata
    report::NamedTuple
end
radial_cutoff(ch::TabulatedRadial)=ch.metadata.cutoff
radial_index_weights(ch::TabulatedRadial)=ch.metadata.index_weights
function radial_ed(ch::TabulatedRadial,r::Real,ps=NamedTuple())
    r>=radial_cutoff(ch) && return (fill(zero(r),length(ch.metadata.spec)),fill(zero(r),length(ch.metadata.spec)))
    first(ch.nodes)<=r<=last(ch.nodes) || throw(DomainError(r,"outside tabulated radial domain"))
    k=clamp(searchsortedlast(ch.nodes,r),1,length(ch.nodes)-1)
    h=ch.nodes[k+1]-ch.nodes[k]; t=(r-ch.nodes[k])/h
    c=ch.coefficients
    v=[((c[4,n,k]*t+c[3,n,k])*t+c[2,n,k])*t+c[1,n,k] for n=1:size(c,2)]
    d=[(3c[4,n,k]*t^2+2c[3,n,k]*t+c[2,n,k])/h for n=1:size(c,2)]
    v,d
end

struct AllocatingJoint{F}
    f::F
end
function (f::AllocatingJoint)(v,d,r,p,c)
    values,derivatives=f.f(r,p,c)
    length(values)==length(v) && length(derivatives)==length(v) || throw(DimensionMismatch("allocating radial output"))
    copyto!(v,values); isnothing(d) || copyto!(d,derivatives)
    nothing
end
struct AllocatingValues{F}
    f::F
end
function (f::AllocatingValues)(v,r,p,c)
    values=f.f(r,p,c)
    length(values)==length(v) || throw(DimensionMismatch("allocating radial output"))
    copyto!(v,values); nothing
end
"""Convenience allocating adapter. Analytic f returns `(values, derivatives)`;
with `derivative=:ad`, f returns values and a scalar ForwardDiff fallback supplies
distance derivatives. Temporary arrays are allocated on each call.
"""
function allocating_radial(f,metadata;config=NamedTuple(),derivative=:analytic)
    CallableRadial(derivative===:ad ? AllocatingValues(f) : AllocatingJoint(f),metadata;config,derivative)
end
function radial_evaluate!(v,d,ch::RadialChannel,r::Real,ps=NamedTuple())
    length(v)==size(ch.coefficients,2) && (isnothing(d) || length(d)==length(v)) || throw(DimensionMismatch("radial output shape"))
    if isnothing(d)
        copyto!(v,radial_values(ch,r))
    else
        values,derivatives=radial_ed(ch,r); copyto!(v,values); copyto!(d,derivatives)
    end
    nothing
end
radial_values(ch::TabulatedRadial,r,ps=NamedTuple())=first(radial_ed(ch,r,ps))
"""Opt-in cubic Hermite snapshot on a raw-distance grid. Check off-grid values
and derivatives; reject tolerances that are not met. C1, generally not C2 at
knots: force-loss/vibration second derivatives have knot discontinuities.
Original custom parameter gradients are not represented by this fixed table.
"""
function tabulate_radial(ch,ps=NamedTuple();domain=(0.001radial_cutoff(ch),radial_cutoff(ch)),nodes=1001,atol=1e-6,derivative_atol=1e-5)
    nodes>=2 && 0<=domain[1]<domain[2]==radial_cutoff(ch) || throw(ArgumentError("tabulation domain must end at the cutoff"))
    xs=collect(range(domain...;length=nodes)); samples=[radial_ed(ch,r,ps) for r in xs]
    maximum(abs,last(samples)[1];init=0.)<=atol && maximum(abs,last(samples)[2];init=0.)<=derivative_atol ||
        throw(ArgumentError("tabulation's zero cutoff extension requires vanishing cutoff values and first derivatives"))
    nr=length(first(samples)[1]); coeff=zeros(4,nr,nodes-1)
    for k=1:nodes-1,n=1:nr
        h=xs[k+1]-xs[k]; y0=samples[k][1][n]; y1=samples[k+1][1][n]
        d0=samples[k][2][n]*h; d1=samples[k+1][2][n]*h
        coeff[:,n,k]=[y0,d0,3(y1-y0)-2d0-d1,2(y0-y1)+d0+d1]
    end
    metadata=ch.metadata
    table=TabulatedRadial(xs,coeff,metadata,NamedTuple())
    ev=0.; ed=0.
    for k=1:nodes-1,fraction in (0.2113248654,0.6180339887)
        r=xs[k]+fraction*(xs[k+1]-xs[k]); v,d=radial_ed(ch,r,ps); vt,dt=radial_ed(table,r)
        ev=max(ev,maximum(abs,v-vt;init=0.)); ed=max(ed,maximum(abs,d-dt;init=0.))
    end
    ev<=atol && ed<=derivative_atol || throw(ArgumentError("tabulation tolerances exceeded: value=$ev derivative=$ed"))
    report=(value_error=ev,derivative_error=ed,atol=Float64(atol),derivative_atol=Float64(derivative_atol),
        precision="Float64",fixed_parameter_snapshot=true,continuity="C1; second derivatives may jump at knots")
    (;channel=TabulatedRadial(xs,coeff,metadata,report),report)
end

const RADIAL_REGISTRY=Dict{String,Any}()
"""Register a versioned, trusted radial factory `(config, metadata) -> evaluator`.
Extensions must register their ID before loading its checkpoints.
"""
function register_radial!(id::AbstractString,factory)
    haskey(RADIAL_REGISTRY,id) && throw(ArgumentError("radial ID already registered: $id"))
    RADIAL_REGISTRY[String(id)]=factory
    nothing
end
register_radial!("PotentialFrameWork.smooth_exponential.v1",(c,m)->SmoothExponential())

function source_radial(source::AbstractString,entrypoint::AbstractString,metadata::RadialMetadata;
                       config=NamedTuple(),dependencies=Dict{String,String}(),allow_code=false,namespace_cache=nothing)
    allow_code || throw(ArgumentError("custom Julia source loading is disabled; pass allow_code=true only for trusted code"))
    for (name,version) in dependencies
        path=Base.find_package(name)
        isnothing(path) && throw(ArgumentError("missing custom radial dependency $name; install it explicitly"))
        project=TOML.parsefile(joinpath(dirname(dirname(path)),"Project.toml"))
        get(project,"version",nothing)==version || throw(ArgumentError("unsupported $name dependency version; expected $version"))
    end
    cache_key=bytes2hex(SHA.sha256(source))
    cached=!isnothing(namespace_cache) && haskey(namespace_cache,cache_key)
    namespace=cached ? namespace_cache[cache_key] : Module(gensym(:ACERadial))
    try
        cached || Base.include_string(namespace,String(source),"saved_ace_radial.jl")
        isnothing(namespace_cache) || (namespace_cache[cache_key]=namespace)
    catch e
        throw(ArgumentError("custom radial source failed: "*sprint(showerror,e)))
    end
    name=Symbol(entrypoint)
    isdefined(namespace,name) || throw(ArgumentError("missing custom radial entrypoint $entrypoint"))
    factory=Base.invokelatest(getfield,namespace,name)
    f=try
        Base.invokelatest(factory,config,metadata)
    catch e
        throw(ArgumentError("custom radial factory failed: "*sprint(showerror,e)))
    end
    record=(route="source",interface_version=1,source=String(source),entrypoint=String(entrypoint),
        source_sha256=bytes2hex(SHA.sha256(source)),dependencies=Dict{String,String}(dependencies))
    CallableRadial(f,metadata;config,persistence=record)
end
has_source(ch)=false
has_source(ch::CallableRadial)=!isnothing(ch.persistence) && get(ch.persistence,:route,"")=="source"
has_source(ch::SubsetRadial)=has_source(ch.parent)
has_source(m::ResolvedModel)=any(has_source,m.radial) || any(has_source,m.pair_radial)

radial_definition_id(ch::CallableRadial)=(ch.metadata,ch.config,isnothing(ch.persistence) ? string(typeof(ch.evaluator)) : ch.persistence)
radial_definition_id(ch::TabulatedRadial)=(ch.metadata,ch.nodes,ch.coefficients)
radial_definition_id(ch::SubsetRadial)=(radial_definition_id(ch.parent),ch.indices)
require_radial_ad(ch::RadialChannel)=nothing
require_radial_ad(ch::TabulatedRadial)=nothing
require_radial_ad(ch::SubsetRadial)=require_radial_ad(ch.parent)
function require_radial_ad(ch::CallableRadial)
    ch.metadata.ad || throw(ArgumentError("custom radial did not declare training AD support; generic second/mixed derivatives are required"))
end

radial_parameter_seed(ps)=parameter_seed(ps)
radial_parameter_seed(::NamedTuple{()})=0

radial_ad_probe(ch,ps,r)=(checked=false,passed=false,reason="load ForwardDiff to check local second/mixed radial derivatives")
