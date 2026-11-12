# Explicit inference conversion. Functional CPU AD never calls this machinery.
convert_data(x::AbstractFloat,T,A)=T(x)
convert_data(x::Rational,T,A)=T(x)
convert_data(x::Number,T,A)=x isa Integer ? deepcopy(x) : throw(ArgumentError("explicit inference conversion does not accept tracked/complex parameters"))
convert_data(x::AbstractArray{<:AbstractFloat},T,A)=A(T.(Array(x)))
convert_data(x::NamedTuple,T,A)=map(v->convert_data(v,T,A),x)
convert_data(x::Tuple,T,A)=map(v->convert_data(v,T,A),x)
convert_data(x::AbstractArray{<:Integer},T,A)=A(copy(Array(x)))
convert_data(x,T,A)=deepcopy(x)
function convert_channel(ch::RadialChannel,::Type{T}) where T
    t=ch.transform
    tr=AgnesiTransform(t.p,t.q,T(t.a),T(t.rin),T(t.r0),T(t.yin),T(t.ycut),T(t.rcut))
    RadialChannel(tr,ch.envelope,T.(ch.envelope_parameters),T.(ch.nodes),T.(ch.coefficients))
end
function convert_channel(ch::CallableRadial,T)
    owned=deepcopy(ch) # copy mutable captures; immutable code/module bindings may be shared
    CallableRadial(owned.evaluator,convert_data(owned.config,T,Array),convert_metadata(owned.metadata,T),owned.persistence)
end
convert_channel(ch::SubsetRadial,T)=SubsetRadial(convert_channel(ch.parent,T),copy(ch.indices))
convert_channel(ch::TabulatedRadial,T)=TabulatedRadial(T.(ch.nodes),T.(ch.coefficients),convert_metadata(ch.metadata,T),deepcopy(ch.report))

"""Host descriptor plus device numerical tables. `source` is an independent CPU
model in selected precision; it is construction metadata, never a GPU fallback.
Use `adapt_model` for complete model/parameter conversion, not field mutation.
"""
struct GPUACEModel{M,D,T}
    source::M
    tables::D
    precision::Type{T}
end

"""Create an independent inference model/parameter bundle. Array selects CPU;
other array constructors select packed storage (e.g. AMDGPU.ROCArray). All
explicit floating parameters are converted together. Float32 → Float64 cannot
recover discarded precision. Custom config types can extend `convert_data`.
"""
function adapt_model(m::ResolvedModel,ps;array_type=Array,precision=eltype(m.reference))
    precision in (Float32,Float64) || throw(ArgumentError("inference precision must be Float32 or Float64"))
    T=precision
    validate_checkpoint_model(m,ps)
    validate_radial_layout(m)
    cpu=ResolvedModel(copy(m.atomic_numbers),T.(m.reference),map(ch->convert_channel(ch,T),m.radial),
        map(ch->convert_channel(ch,T),m.pair_radial),copy(m.radial_spec),copy(m.pair_spec),m.lmax,
        copy(m.density),deepcopy(m.products),SparseMatrixCSC{T,Int}(m.coupling),m.source_sha256,deepcopy(m.policy),m.coordinate_kind,(T===Float32 || m.coefficient_precision===Float32 || parameter_precision(ps)===Float32) ? Float32 : Float64,m.angular_kind)
    validate_radial_layout(cpu)
    pars=convert_data(ps,T,array_type)
    model=array_type===Array ? cpu : GPUACEModel(cpu,pack_model(cpu,array_type),T)
    (;model,parameters=pars)
end
adapt_model(m::GPUACEModel,ps;array_type=Array,precision=m.precision)=
    adapt_model(m.source,map(x->host_data(x),ps);array_type,precision)
host_data(x::AbstractArray)=Array(x)
host_data(x::NamedTuple)=map(host_data,x)
host_data(x)=x

# Ragged spline tables are padded only along model dimensions, never edges.
function pack_radials(channels,A,T)
    if !all(ch->ch isa RadialChannel,channels)
        nr=length(channel_spec(first(channels)))
        return (;kind=Val(:mixed),channels=Tuple(pack_channel(ch,A,T) for ch in channels),nz=size(channels,1),nr)
    end
    nz=size(channels,1); nr=size(first(channels).coefficients,2)
    nk=maximum(ch->length(ch.nodes),channels)
    nodes=zeros(T,nk,nz,nz); coeff=zeros(T,4,nr,nk-1,nz,nz)
    transform=zeros(T,6,nz,nz); powers=zeros(Int,2,nz,nz)
    envelope=zeros(T,5,nz,nz); kinds=zeros(Int,nz,nz); lengths=zeros(Int,nz,nz)
    for iz=1:nz,jz=1:nz
        ch=channels[iz,jz]; t=ch.transform; n=length(ch.nodes)
        nodes[1:n,iz,jz]=ch.nodes; coeff[:,:,1:n-1,iz,jz]=ch.coefficients
        transform[:,iz,jz]=T[t.a,t.rin,t.r0,t.yin,t.ycut,t.rcut]
        powers[:,iz,jz]=[t.p,t.q]; lengths[iz,jz]=n
        ch.envelope in (:PolyEnvelope2sX,:ACE1_PolyEnvelope1sR) || throw(ArgumentError("GPU packing does not support $(ch.envelope)"))
        kinds[iz,jz]=ch.envelope===:PolyEnvelope2sX ? 1 : 2
        envelope[1:length(ch.envelope_parameters),iz,jz]=ch.envelope_parameters
    end
    (;kind=Val(:spline),nodes=A(nodes),coeff=A(coeff),transform=A(transform),powers=A(powers),envelope=A(envelope),kinds=A(kinds),lengths=A(lengths),nr)
end
function pack_model(m,A)
    m.angular_kind===:real_spherical || throw(ArgumentError("GPU ACE does not support solid harmonics"))
    T=eltype(m.reference); nu=maximum(length,m.products;init=0)
    products=zeros(Int,nu,length(m.products)); orders=length.(m.products)
    for j in eachindex(m.products); products[1:orders[j],j]=m.products[j]; end
    C=m.coupling; rows=Int[]; cols=Int[]; vals=T[]
    for j=1:size(C,2), k=nzrange(C,j)
        m.policy.mask[rowvals(C)[k]] || continue
        push!(rows,rowvals(C)[k]); push!(cols,j); push!(vals,nonzeros(C)[k])
    end
    selected=[z in m.policy.species for z in eachindex(m.atomic_numbers)]
    pairmask=[pair_allowed(m,z) for z in eachindex(m.atomic_numbers)]
    range=isnothing(m.policy.orders) ? (0:nu) : m.policy.orders
    (;radial=pack_radials(m.radial,A,T),pair=pack_radials(m.pair_radial,A,T),
      reference=A([reference_value(m,z) for z in eachindex(m.atomic_numbers)]),
      density=A(reshape(collect(Iterators.flatten(m.density)),2,:)),products=A(products),orders=A(orders),
      rows=A(rows),cols=A(cols),coffsets=A(vcat(1,1 .+ cumsum([count(==(j),cols) for j=1:size(C,2)]))),coupling=A(vals),selected=A(selected),pairmask=A(pairmask),
      minorder=first(range),maxorder=last(range),degree=nu,lmax=m.lmax,
      selected_storage=any(selected) && !(first(range)<=0 && last(range)>=nu),
      normalization=A(BlockAngular(m.lmax;precision=T,width=1).normalization))
end

pack_channel(ch::RadialChannel,A,T)=(;kind=Val(:spline_single),data=pack_radials(reshape([ch],1,1),A,T))
function pack_channel(ch::CallableRadial,A,T)
    isbitstype(typeof(ch.evaluator)) || throw(ArgumentError("GPU callable code must be isbits; keep mutable data in explicit config/parameters"))
    ch.evaluator isa ADRadial && throw(ArgumentError("GPU custom radials require explicit first derivatives"))
    (;kind=Val(:callable),evaluator=ch.evaluator,config=convert_data(ch.config,T,A),count=Val(length(ch.metadata.spec)),cutoff=T(radial_cutoff(ch)))
end
pack_channel(ch::SubsetRadial,A,T)=(;kind=Val(:subset),parent=pack_channel(ch.parent,A,T),indices=A(copy(ch.indices)))
pack_channel(ch::TabulatedRadial,A,T)=(;kind=Val(:tabulated),nodes=A(T.(ch.nodes)),coefficients=A(T.(ch.coefficients)),cutoff=T(radial_cutoff(ch)))

convert_metadata(m::RadialMetadata,T)=RadialMetadata(copy(m.spec),T(m.cutoff),deepcopy(m.index_weights),m.simd,m.ad,m.origin,m.envelope,m.channel_dependence)

# Explicit inference-only geometry conversion before topology construction.
struct PrecisionGeometry{S,T} <: AtomsBase.AbstractSystem{3}
    source::S
    positions::Matrix{T}
    cell::Matrix{T}
end
Base.length(s::PrecisionGeometry)=size(s.positions,2)
AtomsBase.position(s::PrecisionGeometry,i::Integer)=SVector{3}(s.positions[:,i])*u"Å"
AtomsBase.cell_vectors(s::PrecisionGeometry)=ntuple(j->SVector{3}(s.cell[:,j])*u"Å",3)
AtomsBase.periodicity(s::PrecisionGeometry)=AtomsBase.periodicity(s.source)
AtomsBase.species(s::PrecisionGeometry,i::Integer)=AtomsBase.species(s.source,i)
function precision_geometry(sys,T)
    n=length(sys); cv=AtomsBase.cell_vectors(sys)
    X=reshape(T[ustrip(u"Å",AtomsBase.position(sys,i)[a]) for a=1:3,i=1:n],3,n)
    H=T[ustrip(u"Å",cv[j][a]) for a=1:3,j=1:3]
    PrecisionGeometry(sys,X,H)
end

# Preserve explicitly parameterized construction of the existing radial wrappers.
CallableRadial{F,C,S}(f,c,m,s) where {F,C,S}=CallableRadial{F,C,S,typeof(m)}(f,c,m,s)
TabulatedRadial{T}(nodes,coefficients,metadata,report) where T=TabulatedRadial{T,typeof(metadata)}(nodes,coefficients,metadata,report)

coherent_float_type(x,T)=true
coherent_float_type(x::AbstractFloat,T)=typeof(x)===T
coherent_float_type(x::AbstractArray{<:AbstractFloat},T)=eltype(x)===T
coherent_float_type(x::NamedTuple,T)=all(v->coherent_float_type(v,T),values(x))
coherent_float_type(x::Tuple,T)=all(v->coherent_float_type(v,T),x)
parameter_precision(x)=Float64
parameter_precision(x::Float32)=Float32
parameter_precision(x::AbstractArray{Float32})=Float32
parameter_precision(x::NamedTuple)=any(v->parameter_precision(v)===Float32,values(x)) ? Float32 : Float64
parameter_precision(x::Tuple)=any(v->parameter_precision(v)===Float32,x) ? Float32 : Float64

function validate_radial_layout(m)
    for (family,spec) in ((m.radial,m.radial_spec),(m.pair_radial,m.pair_spec))
        for ch in family
            if ch isa RadialChannel || ch isa TabulatedRadial
                nodes=ch.nodes; coefficients=ch.coefficients
                length(nodes)>=2 && all(isfinite,nodes) && all(>(0),diff(nodes)) || throw(ArgumentError("invalid or precision-collapsed radial knots"))
                size(coefficients)==(4,length(spec),length(nodes)-1) || throw(DimensionMismatch("radial coefficient/channel shape"))
                all(isfinite,coefficients) || throw(ArgumentError("nonfinite radial coefficients in selected precision"))
            else
                channel_spec(ch)==spec || throw(DimensionMismatch("radial metadata must match the model channel ordering"))
            end
        end
    end
    nothing
end
