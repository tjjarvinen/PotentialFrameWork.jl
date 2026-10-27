# Portable device inference. Vendor packages are deliberately absent here.
import KernelAbstractions as KA
import AcceleratedKernels as AK
using KernelAbstractions: @kernel, @index, @localmem, @synchronize

# AK 0.4 remains supported for the existing SpheriCart 0.2 CPU environment.
# The tested AMDGPU 2.10 environment uses AK 0.5/KA 0.10.
@static if isdefined(AK,:workspace)
    gpu_sort_workspace(a)=AK.workspace(AK.sort!,a;lt=record_less)
    gpu_sort!(a,s)=AK.sort!(a;lt=record_less,workspace=s)
    gpu_reduce_workspace(dst,src)=AK.workspace(AK.mapreducedim!,identity,+,dst,src;init=zero(eltype(src)),overwrite=true)
    gpu_reduce!(dst,src,s)=AK.mapreducedim!(identity,+,dst,src;init=zero(eltype(src)),overwrite=true,workspace=s)
    gpu_scan_workspace(a)=AK.workspace(AK.accumulate!,+,a;init=0)
    gpu_scan!(a,s)=AK.accumulate!(+,a;init=0,workspace=s)
    gpu_max_workspace(a)=AK.workspace(AK.reduce,max,a;init=zero(eltype(a)))
    gpu_max(a,s)=AK.reduce(max,a;init=zero(eltype(a)),workspace=s)
else
    gpu_sort_workspace(a)=similar(a)
    gpu_sort!(a,s)=AK.sort!(a;lt=record_less,temp=s)
    gpu_reduce_workspace(dst,src)=nothing
    gpu_reduce!(dst,src,s)=AK.reduce(+,src;init=zero(eltype(src)),dims=2,temp=dst)
    gpu_scan_workspace(a)=nothing
    gpu_scan!(a,s)=AK.accumulate!(+,a;init=0)
    gpu_max_workspace(a)=nothing
    gpu_max(a,s)=AK.reduce(max,a;init=zero(eltype(a)))
end

function gpu_length(::Type{T},x) where T
    x isa Unitful.Length || throw(ArgumentError("GPU positions/cell require length units"))
    # Nonfinite open cell generators are allowed; positions checked separately.
    value=ustrip(u"Å",x)
    value isa Union{AbstractFloat,Integer,Rational} || throw(ArgumentError("GPU differentiation is unsupported"))
    T(value)
end

"""Host extraction for the GPU single-image boundary. No CPU neighbor building."""
function gpu_inputs(sys,radius,::Type{T},config) where T
    c=T(radius); search=c+T(config.skin)
    isfinite(c) && c>0 && isfinite(search) || throw(ArgumentError("positive finite cutoff required"))
    pb=Tuple(Bool.(AtomsBase.periodicity(sys)))
    cv=AtomsBase.cell_vectors(sys)
    H=[gpu_length(T,cv[j][a]) for a=1:3,j=1:3]
    all(a==j || iszero(H[a,j]) for a=1:3,j=1:3) || throw(ArgumentError("GPU skew cells are unsupported"))
    L=ntuple(a->H[a,a],3)
    all(!pb[a] || (isfinite(L[a]) && L[a]>T(2)*search) for a=1:3) ||
        throw(ArgumentError("GPU requires periodic lengths strictly greater than twice cutoff+skin; multi-image GPU geometry unsupported"))
    n=length(sys)
    config.skin>0 && n>isqrt(typemax(Int)) && throw(OverflowError("GPU candidate offset capacity"))
    all(!pb[a] || L[a]/search<typemax(Int64)÷4 for a=1:3) || throw(OverflowError("GPU periodic cell key overflow"))
    n<=typemax(Int)÷13 || throw(OverflowError("GPU atom offset overflow"))
    x=[gpu_length(T,position(sys,i)[a]) for a=1:3,i=1:n]
    all(isfinite,x) || throw(ArgumentError("GPU requires finite positions"))
    # Leave an integer margin for neighboring keys and Float32 rounding.
    all(abs(v/search)<typemax(Int64)÷4 for v in x) || throw(OverflowError("GPU cell key overflow"))
    counts=ntuple(a->pb[a] ? max(Int64(1),floor(Int64,L[a]/search)) : Int64(1),3)
    labels=[AtomsBase.species(sys,i) for i=1:n]
    ids=reinterpret(UInt64,labels) |> collect
    (;positions=x,identities=ids,species=unique(labels),lengths=L,periodic=pb,counts,radius=c)
end

struct GPUConfig
    strategy::Symbol
    workgroup::Int
    tile::Int
    skin::Float64
    function GPUConfig(;strategy=:central,workgroup=64,tile=4,skin=0)
        strategy in (:central,:candidate) || throw(ArgumentError("unknown GPU strategy"))
        workgroup in (32,64,128,256) || throw(ArgumentError("workgroup must be 32/64/128/256"))
        tile in (1,2,4) || throw(ArgumentError("tile must be 1/2/4"))
        isfinite(skin) && skin>=0 || throw(ArgumentError("skin must be finite and nonnegative"))
        new(strategy,workgroup,tile,skin)
    end
end
struct GPURecord{T}
    key::NTuple{3,Int64}
    xyz::NTuple{3,T}
    species::UInt64
    original::Int
end
record_less(a,b)=isless(a.key,b.key)

mutable struct GPUWorkspace{T,B,R,A,I}
    backend::B
    config::GPUConfig
    records::R
    sorttemp::Any
    x::A
    ends::I
    species::I
    permutation::I
    partials::A
    totals::A
    force::A
    n::Int
    lengths::NTuple{3,T}
    periodic::NTuple{3,Bool}
    counts::NTuple{3,Int64}
    radius::T
    cache::Any
    rebuilds::Int
    metadata::Any
    transfer::Any
    reducetemp::Any
    generation::UInt
end
function GPUWorkspace(backend;precision=Float64,config=GPUConfig())
    precision in (Float32,Float64) || throw(ArgumentError("GPU precision must be Float32 or Float64; GPU AD is unsupported"))
    T=precision
    config.strategy===:candidate && (13*config.tile*config.workgroup+3*config.tile)*sizeof(T)>65536 &&
        throw(ArgumentError("candidate local memory exceeds 64 KiB; reduce tile/workgroup"))
    alloc(S,dims...)=KA.allocate(backend,S,dims...)
    GPUWorkspace(backend,config,alloc(GPURecord{T},0),alloc(GPURecord{T},0),
        alloc(T,3,0),alloc(Int,0),alloc(Int,0),alloc(Int,0),alloc(T,10,0),
        alloc(T,10,1),alloc(T,3,0),0,ntuple(_->one(T),3),(false,false,false),
        (Int64(1),Int64(1),Int64(1)),zero(T),nothing,0,nothing,nothing,nothing,UInt(0))
end

# Cell keys are sparse integer triples, never vacuum-sized dense tables.
@inline function gpu_key(x,L,pb,counts,c)
    @inbounds begin
    ntuple(Val(3)) do a
        @inbounds begin
        if pb[a]
            q=mod(x[a],L[a]); min(unsafe_trunc(Int64,floor(q/(L[a]/counts[a]))),counts[a]-1)
        else
            unsafe_trunc(Int64,floor(x[a]/c))
        end
        end
    end
    end
end
@kernel function gpu_records!(records,x,ids,L,pb,counts,c)
    @inbounds begin
    i=@index(Global)
    xyz=ntuple(a->(@inbounds x[a,i]),Val(3))
    records[i]=GPURecord(gpu_key(xyz,L,pb,counts,c),xyz,ids[i],i)
    end
end
@inline function key_lower(records,key,n)
    @inbounds begin
    lo=1; hi=n+1
    while lo<hi
        mid=(lo+hi)>>>1
        if isless(records[mid].key,key); lo=mid+1; else; hi=mid; end
    end
    lo
    end
end
@inline function key_upper(records,key,n)
    @inbounds begin
    lo=1; hi=n+1
    while lo<hi
        mid=(lo+hi)>>>1
        if !isless(key,records[mid].key); lo=mid+1; else; hi=mid; end
    end
    lo
    end
end
@kernel function gpu_unpack!(x,ends,perm,channels,records,labels,n)
    @inbounds begin
    i=@index(Global)
    rec=records[i]
    for a=1:3; x[a,i]=rec.xyz[a]; end
    perm[i]=rec.original
    ch=0
    for k in eachindex(labels); if labels[k]==rec.species; ch=k; end; end
    channels[i]=ch
    if i==1 || records[i-1].key!=rec.key
        ends[i]=key_upper(records,rec.key,n)-1
    end
    end
end
# Enumerate unique surrounding keys, including the two-cell periodic case.
@inline function adjacent_key(key,t,L,pb,counts)
    @inbounds begin
    off=(mod(t,3)-1,mod(t÷3,3)-1,t÷9-1)
    valid=true
    target=ntuple(Val(3)) do a
        @inbounds begin
        k=key[a]+off[a]
        pb[a] ? mod(k,counts[a]) : k
        end
    end
    for a=1:3
        if pb[a]
            valid &= !(counts[a]==1 && off[a]!=0)
            valid &= !(counts[a]==2 && off[a]==1)
        end
    end
    target,valid
    end
end
@inline gpu_delta(x,i,j,L,pb)=ntuple(a->@inbounds(begin
    d=x[a,j]-x[a,i]
    pb[a] ? d-L[a]*round(d/L[a]) : d
end),Val(3))

# Decode the dependency's verified native eight-byte representation without
# struct reinterpretation, which POCL's compiler cannot lower.
@inline function gpu_species(bits::UInt64)
    @static if Base.ENDIAN_BOM==0x04030201
        z=UInt16(bits & 0xffff)
        neutrons=UInt16((bits>>16) & 0xffff)
        name=UInt32(bits>>32)
    else
        z=UInt16(bits>>48)
        neutrons=UInt16((bits>>32) & 0xffff)
        name=UInt32(bits & 0xffffffff)
    end
    ChemicalSpecies(reinterpret(Int16,z),reinterpret(Int16,neutrons),name)
end

struct GPUPairView{T}
    r::T
    center::ChemicalSpecies
    neighbor::ChemicalSpecies
    ci::Int
    cj::Int
end
NI.distance(b::GPUPairView)=b.r
NI.central_species(b::GPUPairView)=b.center
NI.neighbour_species(b::GPUPairView)=(b.neighbor,)
NI.neighbour_species(b::GPUPairView,k)=b.neighbor
NI.neighbour_count(::GPUPairView)=1
NI.active(::GPUPairView)=true
NI.structural(::GPUPairView)=true
mutable struct GPUPairOutput{T}
    energy::T
    derivative::T
end
NI.set_energy!(o::GPUPairOutput,e)=(o.energy=e; nothing)
NI.set_radial_derivative!(o::GPUPairOutput,d)=(o.derivative=d; nothing)
NI.output_energy(o::GPUPairOutput)=o.energy
NI.output_radial_derivative(o::GPUPairOutput)=o.derivative

@inline function gpu_pair(x,records,channels,i,j,L,pb,c,p,ps,props,center=ntuple(a->(@inbounds x[a,i]),Val(3)),neighbor=ntuple(a->(@inbounds x[a,j]),Val(3)))
    @inbounds begin
    T=eltype(x); z=zero(T)
    d=ntuple(a->@inbounds(begin
        v=neighbor[a]-center[a]; pb[a] ? v-L[a]*round(v/L[a]) : v
    end),Val(3)); r2=sum(abs2,d)
    if i==j || r2>=c*c
        return ntuple(_->z,Val(13))
    end
    r=sqrt(r2)
    b=GPUPairView(r,gpu_species(records[i].species),
        gpu_species(records[j].species),channels[i],channels[j])
    out=GPUPairOutput(z,z)
    NI.evaluate_pair!(out,p,b,ps,nothing,props)
    # Error sentinel is consumed on host only when explicitly checking results.
    bad=(!isfinite(out.energy) || !isfinite(out.derivative) || (r==z && out.derivative!=z))
    coefficient=out.derivative/ifelse(r==z,one(T),r)
    g=ntuple(a->(@inbounds bad ? T(NaN) : coefficient*d[a]),Val(3))
    h=inv(T(2))
    ntuple(Val(13)) do a
        @inbounds begin
        a==1 ? (bad ? T(NaN) : h*out.energy) : a<=4 ? g[a-1] :
            -h*d[mod(a-5,3)+1]*g[(a-5)÷3+1]
        end
    end
    end
end
@inline add13(a::NTuple{13,T},b::NTuple{13,T}) where T = ntuple(k->(@inbounds a[k]+b[k]),Val(13))
@kernel function gpu_central!(force,partials,x,records,ends,channels,perm,L,pb,counts,c,p,ps,props,n)
    @inbounds begin
    i=@index(Global)
    T=eltype(x); acc=ntuple(_->zero(T),Val(13))
    for t=0:26
        key,valid=adjacent_key(records[i].key,t,L,pb,counts)
        if valid
            lo=key_lower(records,key,n)
            if lo<=n && records[lo].key==key
                for j=lo:ends[lo]
                    v=gpu_pair(x,records,channels,i,j,L,pb,c,p,ps,props)
                    acc=add13(acc,v)
                end
            end
        end
    end
    for a=1:3; force[a,perm[i]]=acc[a+1]; end
    partials[1,i]=acc[1]
    for a=1:9; partials[a+1,i]=acc[a+4]; end
    end
end

# One workgroup owns a bounded sorted central tile. All barriers are uniform.
# Candidate lanes load contiguous cell runs and visit every center in that cell
# tile; no center is ever split across workgroups and no force atomics are used.
@kernel function gpu_candidate!(force,partials,x,records,ends,channels,perm,L,pb,counts,c,p,ps,props,n,::Val{B},::Val{C}) where {B,C}
    @inbounds begin
    lane=@index(Local,Linear)
    group=@index(Group,Linear)
    sums=@localmem eltype(x) (13,C,B)
    centers=@localmem eltype(x) (3,C)
    for k=1:C
        i=(group-1)*C+k
        if lane<=3; centers[lane,k]=i<=n ? x[lane,i] : zero(eltype(x)); end
        for a=1:13; sums[a,k,lane]=zero(eltype(x)); end
    end
    @synchronize
    for anchor=1:C
        first=(group-1)*C+anchor
        if first<=n
            key0=records[first].key
            unique_cell=anchor==1 || records[first-1].key!=key0
            if unique_cell
                for t=0:26
                    key,valid=adjacent_key(key0,t,L,pb,counts)
                    lo=key_lower(records,key,n)
                    if valid && lo<=n && records[lo].key==key
                        for j=lo+lane-1:B:ends[lo]
                            candidate=ntuple(a->(@inbounds x[a,j]),Val(3))
                            for k=1:C
                                i=(group-1)*C+k
                                if i<=n && records[i].key==key0
                                    v=gpu_pair(x,records,channels,i,j,L,pb,c,p,ps,props,ntuple(a->(@inbounds centers[a,k]),Val(3)),candidate)
                                    for a=1:13; sums[a,k,lane]+=v[a]; end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    @synchronize
    for level=1:trailing_zeros(B)
        step=B>>level
        if lane<=step
            for k=1:C,a=1:13; sums[a,k,lane]+=sums[a,k,lane+step]; end
        end
        @synchronize
    end
    if lane==1
        for k=1:C
            i=(group-1)*C+k
            if i<=n
                for a=1:3; force[a,perm[i]]=sums[a+1,k,1]; end
                partials[1,i]=sums[1,k,1]
                for a=1:9; partials[a+1,i]=sums[a+4,k,1]; end
            end
        end
    end
    end
end

function gpu_resize!(w::GPUWorkspace{T},n) where T
    n==w.n && return w
    alloc(S,dims...)=KA.allocate(w.backend,S,dims...)
    w.records=alloc(GPURecord{T},n); w.sorttemp=gpu_sort_workspace(w.records)
    w.x=alloc(T,3,n); w.ends=alloc(Int,n); w.species=alloc(Int,n)
    w.permutation=alloc(Int,n); w.force=alloc(T,3,n); w.partials=alloc(T,10,n)
    w.reducetemp=gpu_reduce_workspace(w.totals,w.partials)
    w.n=n
    w.cache=nothing
    w
end

"""Rebuild sparse sorted cells from original-order device numeric coordinates.
This mutating inference path has no GPU AD support. Parameters are supplied anew
on every execution. Returned arrays own storage, independent of the workspace.
"""
function gpu_execute!(w::GPUWorkspace{T},positions,identities,labels,p,ps;
        properties=(AC.Energy(),AC.Forces(),AC.Virial())) where T
    props=normalize_properties(properties)
    n=size(positions,2)
    size(positions)==(3,n) && length(identities)==n || throw(DimensionMismatch())
    eltype(positions)===T || throw(ArgumentError("positions must have selected GPU precision; GPU AD unsupported"))
    gpu_resize!(w,n)
    if n>0 && !isempty(props)
        if w.config.skin>0
            gpu_skin!(w,positions,identities,labels)
            args=(w.force,w.partials,w.x,w.records,w.species,w.permutation,w.lengths,w.periodic,w.radius,p,ps,props,w.cache.offsets,w.cache.neighbors)
            if w.config.strategy===:central
                gpu_buffered_central!(w.backend,w.config.workgroup)(args...;ndrange=n)
            else
                B=w.config.workgroup; C=w.config.tile
                gpu_buffered_candidate!(w.backend,B)(args...,n,Val(B),Val(C);ndrange=Base.checked_mul(cld(n,C),B))
            end
        else
            gpu_build_cells!(w,positions,identities,labels)
            args=(w.force,w.partials,w.x,w.records,w.ends,w.species,w.permutation,w.lengths,w.periodic,w.counts,w.radius,p,ps,props,n)
            if w.config.strategy===:central
                gpu_central!(w.backend,w.config.workgroup)(args...;ndrange=n)
            else
                B=w.config.workgroup; C=w.config.tile
                gpu_candidate!(w.backend,B)(args...,Val(B),Val(C);ndrange=Base.checked_mul(cld(n,C),B))
            end
        end
        gpu_reduce!(w.totals,w.partials,w.reducetemp)
        # Pair callbacks signal an invalid/singular origin through the energy
        # partial even for force-only requests. Read one scalar, never forces,
        # so the numerical entrypoint rejects invalid mathematics explicitly.
        status=only(Array(view(w.totals,1:1,1:1)))
        isfinite(status) || throw(DomainError(status,"singular origin or nonfinite GPU radial callback"))
    else
        fill!(w.totals,zero(T)); fill!(w.force,zero(T))
    end
    merge(wants(AC.Energy,props) ? (;energy=copy(view(w.totals,1:1,:))) : (;),
        wants(AC.Forces,props) ? (;forces=copy(w.force)) : (;),
        wants(AC.Virial,props) ? (;virial=reshape(copy(view(w.totals,2:10,:)),3,3)) : (;))
end

@kernel function gpu_displacements!(distance2,positions,reference,identities,oldids,labels,oldlabels)
    @inbounds begin
    i=@index(Global)
    d=sum(a->(@inbounds abs2(positions[a,i]-reference[a,i])),1:3)
    identities[i]==oldids[i] || (d=eltype(positions)(Inf))
    if i==1
        for k in eachindex(labels)
            labels[k]==oldlabels[k] || (d=eltype(positions)(Inf))
        end
    end
    distance2[i]=d
    end
end
@kernel function gpu_update_positions!(x,positions,permutation)
    @inbounds begin
    i=@index(Global)
    for a=1:3; x[a,i]=positions[a,permutation[i]]; end
    end
end
# Count first, prefix-scan, then fill exact conservative CSR storage. No truncation.
@kernel function gpu_connectivity!(offsets,neighbors,x,records,ends,L,pb,counts,search,n,::Val{Fill}) where Fill
    @inbounds begin
    i=@index(Global)
    count=0
    for t=0:26
        key,valid=adjacent_key(records[i].key,t,L,pb,counts)
        if valid
            lo=key_lower(records,key,n)
            if lo<=n && records[lo].key==key
                for j=lo:ends[lo]
                    d=gpu_delta(x,i,j,L,pb)
                    if i!=j && sum(abs2,d)<=search*search
                        count+=1
                        if Fill; neighbors[(i==1 ? 0 : offsets[i-1])+count]=j; end
                    end
                end
            end
        end
    end
    if !Fill; offsets[i]=count; end
    end
end
@kernel function gpu_buffered_central!(force,partials,x,records,channels,perm,L,pb,c,p,ps,props,offsets,neighbors)
    @inbounds begin
    i=@index(Global)
    acc=ntuple(_->zero(eltype(x)),Val(13))
    lo=(i==1 ? 0 : offsets[i-1])+1
    for k=lo:offsets[i]
        v=gpu_pair(x,records,channels,i,neighbors[k],L,pb,c,p,ps,props)
        acc=add13(acc,v)
    end
    for a=1:3; force[a,perm[i]]=acc[a+1]; end
    partials[1,i]=acc[1]
    for a=1:9; partials[a+1,i]=acc[a+4]; end
    end
end
@kernel function gpu_buffered_candidate!(force,partials,x,records,channels,perm,L,pb,c,p,ps,props,offsets,neighbors,n,::Val{B},::Val{C}) where {B,C}
    @inbounds begin
    lane=@index(Local,Linear)
    group=@index(Group,Linear)
    sums=@localmem eltype(x) (13,C,B)
    centers=@localmem eltype(x) (3,C)
    for k=1:C
        i=(group-1)*C+k
        if lane<=3; centers[lane,k]=i<=n ? x[lane,i] : zero(eltype(x)); end
        for a=1:13; sums[a,k,lane]=zero(eltype(x)); end
    end
    @synchronize
    for k=1:C
        i=(group-1)*C+k
        if i<=n
            lo=(i==1 ? 0 : offsets[i-1])+1
            for entry=lo+lane-1:B:offsets[i]
                v=gpu_pair(x,records,channels,i,neighbors[entry],L,pb,c,p,ps,props,ntuple(a->(@inbounds centers[a,k]),Val(3)))
                for a=1:13; sums[a,k,lane]+=v[a]; end
            end
        end
    end
    @synchronize
    for level=1:trailing_zeros(B)
        step=B>>level
        if lane<=step
            for k=1:C,a=1:13; sums[a,k,lane]+=sums[a,k,lane+step]; end
        end
        @synchronize
    end
    if lane==1
        for k=1:C
            i=(group-1)*C+k
            if i<=n
                for a=1:3; force[a,perm[i]]=sums[a+1,k,1]; end
                partials[1,i]=sums[1,k,1]
                for a=1:9; partials[a+1,i]=sums[a+4,k,1]; end
            end
        end
    end
    end
end

function gpu_build_cells!(w,positions,identities,labels)
    search=w.radius+typeof(w.radius)(w.config.skin)
    gpu_records!(w.backend,w.config.workgroup)(w.records,positions,identities,w.lengths,w.periodic,w.counts,search;ndrange=w.n)
    gpu_sort!(w.records,w.sorttemp)
    gpu_unpack!(w.backend,w.config.workgroup)(w.x,w.ends,w.permutation,w.species,w.records,labels,w.n;ndrange=w.n)
    w.rebuilds+=1
end
function gpu_skin!(w::GPUWorkspace{T},positions,identities,labels) where T
    # The numerical caller must invalidate cache on metadata/species changes.
    # Raw-coordinate endpoint displacement bounds torus displacement. Wrapped
    # crossings rebuild conservatively; no winding is inferred or required.
    signature=(w.lengths,w.periodic,w.counts,w.radius,w.config.skin,length(labels),w.n)
    rebuild=isnothing(w.cache) || w.cache.signature!=signature
    if !rebuild
        gpu_displacements!(w.backend,w.config.workgroup)(w.cache.distance2,positions,w.cache.reference,identities,w.cache.identities,labels,w.cache.labels;ndrange=w.n)
        max2=gpu_max(w.cache.distance2,w.cache.maxworkspace) # identified host scalar read
        rebuild=T(2)*sqrt(max2)>=T(w.config.skin)*(one(T)-T(16)*eps(T))
    end
    if rebuild
        gpu_build_cells!(w,positions,identities,labels)
        old=w.cache
        reusable=!isnothing(old) && length(old.offsets)==w.n
        offsets=reusable ? old.offsets : KA.allocate(w.backend,Int,w.n)
        empty=KA.allocate(w.backend,Int,0)
        search=w.radius+T(w.config.skin)
        # Inflate radius slightly to conservatively include rounding at the bound.
        search=nextfloat(search,16)
        gpu_connectivity!(w.backend,w.config.workgroup)(offsets,empty,w.x,w.records,w.ends,w.lengths,w.periodic,w.counts,search,w.n,Val(false);ndrange=w.n)
        scanworkspace=reusable ? old.scanworkspace : gpu_scan_workspace(offsets)
        gpu_scan!(offsets,scanworkspace)
        count=only(Array(view(offsets,w.n:w.n))) # exact capacity read on rebuild
        count>=0 || throw(OverflowError("GPU connectivity offset overflow"))
        neighbors=reusable && length(old.neighbors)>=count ? old.neighbors : KA.allocate(w.backend,Int,count)
        gpu_connectivity!(w.backend,w.config.workgroup)(offsets,neighbors,w.x,w.records,w.ends,w.lengths,w.periodic,w.counts,search,w.n,Val(true);ndrange=w.n)
        reference=reusable ? old.reference : similar(positions)
        distance2=reusable ? old.distance2 : KA.allocate(w.backend,T,w.n)
        copyto!(reference,positions)
        maxworkspace=reusable ? old.maxworkspace : gpu_max_workspace(distance2)
        w.cache=(;offsets,neighbors,signature,scanworkspace,maxworkspace,identities=copy(identities),labels=copy(labels),reference,distance2)
    else
        gpu_update_positions!(w.backend,w.config.workgroup)(w.x,positions,w.permutation;ndrange=w.n)
    end
    nothing
end
