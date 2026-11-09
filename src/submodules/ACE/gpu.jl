import KernelAbstractions as KA
using KernelAbstractions: @kernel, @index, @atomic, @synchronize, @localmem, @private

# One scalar channel, exact imported cubic + transform + envelope and d/dr.
@inline gpu_radial(ch,r,n,iz,jz,ps=NamedTuple())=gpu_radial(ch.kind,ch,r,n,iz,jz,ps)
@inline function gpu_radial(::Val{:spline},ch,r,n,iz,jz,ps)
    @inbounds begin
        a,rin,r0,yin,ycut,rc=ntuple(k->ch.transform[k,iz,jz],6)
        r>=rc && return (zero(r),zero(r))
        p=ch.powers[1,iz,jz]; q=ch.powers[2,iz,jz]; h=q-p
        s=(r-rin)/(r0-rin)
        y=r<=rin ? one(r) : inv(one(r)+a*s^q/(one(r)+s^h))
        x=clamp(-one(r)+2*(y-yin)/(ycut-yin),-one(r),one(r))
        dx=zero(r)
        if r>rin && -one(r)<x<one(r)
            den=one(r)+s^h; u=s^q/den
            du=(q*s^(q-1)/den-(h==0 ? zero(r) : h*s^(q+h-1)/den^2))/(r0-rin)
            dx=-2*a*du/((one(r)+a*u)^2*(ycut-yin))
        end
        ep=ntuple(k->ch.envelope[k,iz,jz],5)
        e=zero(r); de=zero(r)
        if ch.kinds[iz,jz]==1
            if ep[1]<x<ep[2]
                pe=Int(ep[3]); qe=Int(ep[4]); left=x-ep[1]; right=ep[2]-x
                e=ep[5]*left^pe*right^qe
                de=ep[5]*(pe*left^(pe-1)*right^qe-qe*left^pe*right^(qe-1))*dx
            end
        else
            pe=Int(ep[3]); sr=r/ep[2]; sc=ep[1]/ep[2]
            e=sr^(-pe)-sc^(-pe)+pe*sc^(-pe-1)*(sr-sc)
            de=pe/ep[2]*(sc^(-pe-1)-sr^(-pe-1))
        end
        lo=1; hi=ch.lengths[iz,jz]
        while lo<hi
            mid=(lo+hi+1)>>>1
            if ch.nodes[mid,iz,jz]<=x; lo=mid; else; hi=mid-1; end
        end
        k=min(lo,ch.lengths[iz,jz]-1)
        step=ch.nodes[k+1,iz,jz]-ch.nodes[k,iz,jz]; t=(x-ch.nodes[k,iz,jz])/step
        c=ntuple(j->ch.coeff[j,n,k,iz,jz],4)
        v=((c[4]*t+c[3])*t+c[2])*t+c[1]
        dv=(3*c[4]*t*t+2*c[3]*t+c[2])/step
        (e*v,de*v+e*dv*dx)
    end
end

# Same real-harmonic recurrence as angular_blocks.jl, one requested component.
# Constant-size Cartesian jets avoid device allocation and polar singularities.
@inline function gpu_angular(R,midx,F)
    T=typeof(R[1]); r=sqrt(sum(abs2,R)); u=ntuple(k->R[k]/r,3)
    jets=ntuple(k->AngularJet(u[k],ntuple(a->((a==k ? one(T) : zero(T))-u[k]*u[a])/r,3)),3)
    x,y,z=jets; r2=x*x+y*y+z*z
    l=isqrt(midx-1); m=midx-(l*l+l+1); ma=abs(m)
    jone=AngularJet(one(T),(zero(T),zero(T),zero(T)))
    jzero=AngularJet(zero(T),(zero(T),zero(T),zero(T)))
    cs=jone; sn=jzero
    for k=1:ma
        old=cs; cs=cs*x-sn*y; sn=sn*x+old*y
    end
    ma==0 && (cs=cs/sqrt(T(2)))
    qmm=jone
    for k=1:ma; qmm=-(2*k-1)*qmm; end
    q=qmm
    if l>ma
        prev=qmm; q=(2*ma+1)*(z*qmm)
        for k=ma+2:l
            nxt=((2*k-1)*(z*q)-(k+ma-1)*(r2*prev))/(k-ma)
            prev=q; q=nxt
        end
    end
    @inbounds F[l+1,ma+1]*(q*(m<0 ? sn : cs))
end

@inline function gpu_scale(ps,key,n,iz,jz,r)
    hasproperty(ps,key) ? (@inbounds getproperty(ps,key).scale[n,iz,jz]) : one(r)
end

# Shared edge mathematics for both density mappings. The stride selects bounded
# neighbor lanes; inactive channel/tail lanes never normalize an invalid direction.
@inline function gpu_density_partial(x,channels,L,pb,c,offsets,neighbors,tables,ps,i,iz,a,firstedge,stride)
    @inbounds begin
        n=tables.density[1,a]; m=tables.density[2,a]
        u=zero(eltype(x)); v=u
        for k=firstedge:stride:offsets[i]
            j=neighbors[k]; R=NL.gpu_delta(x,i,j,L,pb); r=sqrt(sum(abs2,R))
            if r<c
                if iszero(r); u=eltype(x)(NaN); v=u
                else
                    jz=channels[j]
                    rv,_=gpu_radial(tables.radial,r,n,iz,jz,radial_parameters(ps,:custom_radial))
                    value=rv*gpu_scale(ps,:radial,n,iz,jz,r)*gpu_angular(R,m,tables.normalization).v
                    if tables.selected_storage && tables.selected[jz]; v+=value; else; u+=value; end
                end
            end
        end
        (u,v)
    end
end
@kernel function ace_density_baseline!(U,V,x,channels,perm,L,pb,c,offsets,neighbors,tables,ps,centers,firstsite,count)
    @inbounds begin
        slot=@index(Global,Linear)
        if slot<=count
            i=centers[firstsite+slot-1]; iz=channels[i]
            firstedge=(i==1 ? 0 : offsets[i-1])+1
            for a=1:size(U,1)
                u,v=gpu_density_partial(x,channels,L,pb,c,offsets,neighbors,tables,ps,i,iz,a,firstedge,1)
                U[a,slot]=u
                if size(V,1)>0; V[a,slot]=v; end
            end
        end
    end
end

# One 64-thread workgroup/site, eight channels by eight neighbor lanes. Only this
# bounded tile occupies shared memory, regardless of basis/neighborhood size.
# Each channel's lane zero owns its density; barriers complete cooperative sums.
@kernel function ace_density!(U,V,x,channels,perm,L,pb,c,offsets,neighbors,tables,ps,centers,firstsite,count)
    state=@private Int (5,) # site, center, species, channel lane, neighbor lane
    slot=@index(Group,Linear); lane=@index(Local,Linear)
    state[1]=slot; state[2]=centers[firstsite+slot-1]; state[3]=channels[state[2]]
    state[4]=(lane-1)%8+1; state[5]=(lane-1)÷8+1
    partial=@localmem eltype(x) (2,8,8)
    @inbounds for firstchannel=1:8:size(U,1)
        a=firstchannel+state[4]-1
        u=zero(eltype(x)); v=u
        if a<=size(U,1)
            firstedge=(state[2]==1 ? 0 : offsets[state[2]-1])+state[5]
            u,v=gpu_density_partial(x,channels,L,pb,c,offsets,neighbors,tables,ps,state[2],state[3],a,firstedge,8)
        end
        partial[1,state[4],state[5]]=u; partial[2,state[4],state[5]]=v
        @synchronize
        for step in (4,2,1)
            if state[5]<=step
                partial[1,state[4],state[5]]+=partial[1,state[4],state[5]+step]
                partial[2,state[4],state[5]]+=partial[2,state[4],state[5]+step]
            end
            @synchronize
        end
        a=firstchannel+state[4]-1
        if state[5]==1 && a<=size(U,1)
            U[a,state[1]]=partial[1,state[4],1]
            if size(V,1)>0; V[a,state[1]]=partial[2,state[4],1]; end
        end
        @synchronize
    end
end

# Polynomial coefficients in U+tV, truncated to the largest retained order.
# Fixed compile-time degree; no division by densities, repeated/zero factors safe.
@inline function gpu_product(U,V,slot,products,j,nu,skip,::Val{D},lo,hi) where D
    T=eltype(U)
    if size(V,1)==0
        lo<=0<=hi || return zero(T)
        result=one(T)
        @inbounds for t=1:nu
            if t!=skip; result*=U[products[t,j],slot]; end
        end
        return result
    end
    poly=MVector{D+1,T}(ntuple(_->zero(T),D+1)); poly[1]=one(T)
    @inbounds for t=1:nu
        t==skip && continue
        a=products[t,j]; u=U[a,slot]; v=V[a,slot]
        for k=D:-1:0
            poly[k+1]=u*poly[k+1]+(k>0 ? v*poly[k] : zero(T))
        end
    end
    value=zero(T)
    @inbounds for k=max(0,lo):min(D,hi); value+=poly[k+1]; end
    value
end
@kernel function ace_weights!(weights,tables,ps)
    @inbounds begin
        index=@index(Global,Linear); np=size(weights,1); j=(index-1)%np+1; iz=(index-1)÷np+1
        w=zero(eltype(weights))
        for k=tables.coffsets[j]:tables.coffsets[j+1]-1
            w+=tables.coupling[k]*ps.WB[tables.rows[k],iz]
        end
        weights[j,iz]=w
    end
end

# Site-exclusive product reverse avoids product-adjoint races. Density and reverse
# use separate global arrays. A product-parallel alternative can share these tables.
@inline function gpu_site_products!(slot,partials,dU,dV,U,V,weights,channels,tables,centers,firstsite,count,derivatives,energy,::Val{D}) where D
    @inbounds begin
            i=centers[firstsite+slot-1]; iz=channels[i]
            e=energy ? tables.reference[iz] : zero(eltype(U))
            if derivatives; for a=1:size(U,1); dU[a,slot]=zero(eltype(U)); if size(dV,1)>0; dV[a,slot]=zero(eltype(U)); end; end; end
            for j in eachindex(tables.orders)
                nu=tables.orders[j]; w=weights[j,iz]
                if energy; e+=w*gpu_product(U,V,slot,tables.products,j,nu,0,Val(D),tables.minorder,tables.maxorder); end
                if derivatives
                    for t=1:nu
                        a=tables.products[t,j]
                        dU[a,slot]+=w*gpu_product(U,V,slot,tables.products,j,nu,t,Val(D),tables.minorder,tables.maxorder)
                        if size(dV,1)>0; dV[a,slot]+=w*gpu_product(U,V,slot,tables.products,j,nu,t,Val(D),tables.minorder-1,tables.maxorder-1); end
                    end
                end
            end
            partials[1,i]=e
    end
end
@kernel function ace_products!(partials,dU,dV,U,V,weights,channels,tables,centers,firstsite,count,derivatives,energy,::Val{D}) where D
    slot=@index(Global,Linear)
    if slot<=count
        gpu_site_products!(slot,partials,dU,dV,U,V,weights,channels,tables,centers,firstsite,count,derivatives,energy,Val(D))
    end
end

# Site-exclusive edge traversal. Neighbor force scatter uses atomics; central
# contribution is reduced in registers before its atomic write. Force and virial
# come from each site's own adjoints, never an assumed opposite directed edge.
@inline function gpu_edge_partial(force,dU,dV,x,channels,perm,L,pb,c,offsets,neighbors,tables,ps,i,iz,slot,firstedge,stride,energy,forces,virial)
    T=eltype(x)
    @inbounds begin
            fc=zero(SVector{3,T}); W=zero(SMatrix{3,3,T,9}); ep=zero(T)
            for k=firstedge:stride:offsets[i]
                j=neighbors[k]; jz=channels[j]; R=NL.gpu_delta(x,i,j,L,pb); r=sqrt(sum(abs2,R))
                if r<c
                    if iszero(r)
                        ep=T(NaN)
                    else
                        g=zero(SVector{3,T})
                        if forces || virial
                            for a=1:size(dU,1)
                                n=tables.density[1,a]; m=tables.density[2,a]
                                rv,rd= gpu_radial(tables.radial,r,n,iz,jz,radial_parameters(ps,:custom_radial))
                                scale=gpu_scale(ps,:radial,n,iz,jz,r)
                                Y=gpu_angular(R,m,tables.normalization)
                                adj=size(dV,1)>0 && tables.selected[jz] ? dV[a,slot] : dU[a,slot]
                                g+=adj*scale*(rd*Y.v*SVector(R)/r+rv*SVector(Y.g))
                            end
                        end
                        if tables.pairmask[jz]
                            for n=1:tables.pair.nr
                                q,dq=gpu_radial(tables.pair,r,n,iz,jz,radial_parameters(ps,:custom_pair))
                                w=ps.Wpair[n,iz]*gpu_scale(ps,:pair,n,iz,jz,r)
                                if energy; ep+=w*q; end
                                if forces || virial; g+=w*dq*SVector(R)/r; end
                            end
                        end
                        if (forces || virial) && !all(isfinite,g); ep=T(NaN); end
                        if forces
                            fc+=g
                            for a=1:3; @atomic force[a,perm[j]]+=-g[a]; end
                        end
                        if virial; W-=SVector(R)*transpose(g); end
                    end
                end
            end
        (ep,fc[1],fc[2],fc[3],Tuple(W)...)
    end
end
@inline function gpu_site_edges!(slot,force,partials,dU,dV,x,channels,perm,L,pb,c,offsets,neighbors,tables,ps,centers,firstsite,count,energy,forces,virial)
    @inbounds begin
        i=centers[firstsite+slot-1]; iz=channels[i]
        firstedge=(i==1 ? 0 : offsets[i-1])+1
        v=gpu_edge_partial(force,dU,dV,x,channels,perm,L,pb,c,offsets,neighbors,tables,ps,i,iz,slot,firstedge,1,energy,forces,virial)
        partials[1,i]+=v[1]
        if forces; for a=1:3; @atomic force[a,perm[i]]+=v[a+1]; end; end
        if virial; for a=1:9; partials[a+1,i]=v[a+4]; end; end
    end
end
@kernel function ace_edges!(force,partials,dU,dV,x,channels,perm,L,pb,c,offsets,neighbors,tables,ps,centers,firstsite,count,energy,forces,virial)
    slot=@index(Global,Linear)
    if slot<=count
        gpu_site_edges!(slot,force,partials,dU,dV,x,channels,perm,L,pb,c,offsets,neighbors,tables,ps,centers,firstsite,count,energy,forces,virial)
    end
end

mutable struct GPUExecutor{M,P,W,A,I,C}
    model::M
    parameters::P
    workspace::W
    positions::A
    identities::I
    labels::I
    centers::C
    sortedcenters::C
    inverse::C
    U::A
    V::A
    dU::A
    dV::A
    weights::A
    batch_size::Int
    strategy::Symbol
    properties::Tuple
end

automatic_gpu_strategy(model,reference,n)=:staged

"""Prepare resident scalar ACE inference. Explicit model adaptation must precede
this call. `:staged` assigns density components to threads; `:baseline` assigns
one thread/site per stage. `:auto` uses measured backend preparation rules,
otherwise staged; see docs/gpu_ace.md.
Geometry is orthorhombic single-image. Reprepare after cell/species changes.
Outputs of `execute_gpu!` own device storage; `gpu_result` returns owning host data.
"""
function GPUExecutor(sys,m::GPUACEModel,ps;batch_size=max(1,length(sys)),strategy=:auto,config=NL.GPUConfig(),centers=nothing)
    strategy in (:auto,:staged,:baseline,:compact) || throw(ArgumentError("GPU ACE strategy must be :auto, :staged, :baseline or :compact"))
    batch_size>0 || throw(ArgumentError("positive site batch size required"))
    T=m.precision; backend=KA.get_backend(m.tables.reference)
    check_execution_device(m.tables.reference)
    check_placement(m.tables,m.tables.reference)
    check_placement(ps,m.tables.reference)
    nz=length(m.source.atomic_numbers)
    size(ps.WB)==(size(m.source.coupling,1),nz) || throw(DimensionMismatch("GPU ACE WB shape"))
    size(ps.Wpair)==(length(m.source.pair_spec),nz) || throw(DimensionMismatch("GPU ACE Wpair shape"))
    for (key,nr) in ((:radial,length(m.source.radial_spec)),(:pair,length(m.source.pair_spec)))
        hasproperty(ps,key) && size(getproperty(ps,key).scale)!=(nr,nz,nz) && throw(DimensionMismatch("GPU ACE radial scale shape"))
    end
    check_gpu_data(ps,backend,T)
    input=NL.gpu_inputs(sys,NI.search_cutoff(m.source),T,config)
    for s in input.species; element_channel(m.source,s); end
    w=NL.GPUWorkspace(backend;precision=T,config)
    w.lengths=input.lengths; w.periodic=input.periodic; w.counts=input.counts; w.radius=input.radius
    NL.gpu_resize!(w,length(sys))
    upload(a)=begin buffer=KA.allocate(backend,eltype(a),size(a)...); copyto!(buffer,a); buffer end
    labels=collect(reinterpret(UInt64,ChemicalSpecies.(m.source.atomic_numbers)))
    n=length(sys); b=min(batch_size,max(n,1)); nd=length(m.source.density)
    alloc(dims...)=KA.allocate(backend,T,dims...)
    selected=collect(center_indices(centers,n))
    GPUExecutor(m,deepcopy(ps),w,upload(input.positions),upload(input.identities),upload(labels),upload(selected),upload(selected),upload(zeros(Int,n)),
        alloc(nd,b),alloc(m.tables.selected_storage ? nd : 0,b),alloc(nd,b),alloc(m.tables.selected_storage ? nd : 0,b),alloc(length(m.source.products),length(labels)),b,
        strategy===:auto ? automatic_gpu_strategy(m,m.tables.reference,n) : strategy,())
end
function execute_gpu!(ex::GPUExecutor;parameters=nothing,properties=(AC.Energy(),AC.Forces(),AC.Virial()),owning=true,diagnostics=nothing)
    props=NL.normalize_properties(properties); w=ex.workspace; m=ex.model; t=m.tables; T=m.precision
    check_execution_device(m.tables.reference)
    size(ex.positions)==(3,w.n) && length(ex.identities)==w.n || throw(DimensionMismatch("use prepare_executor! after input capacity changes"))
    check_placement(ex.positions,m.tables.reference)
    if !isnothing(parameters)
        # Explicit parameter update boundary: no precontracted stale readout.
        check_gpu_data(parameters,w.backend,T)
        check_placement(parameters,m.tables.reference)
        check_update_shapes(ex.parameters,parameters)
        ex.parameters=update_gpu_data!(ex.parameters,parameters)
    end
    energy=NL.wants(AC.Energy,props); forces=NL.wants(AC.Forces,props); virial=NL.wants(AC.Virial,props)
    ex.properties=props
    fill!(w.force,zero(T)); fill!(w.partials,zero(T)); fill!(w.totals,zero(T))
    stamp=isnothing(diagnostics) ? UInt64(0) : (KA.synchronize(w.backend); time_ns())
    n=w.n
    if n>0 && !isempty(props)
        # skin=0 forces fresh cells/connectivity; positive skin reuses the same
        # directed CSR machinery with the existing displacement validity check.
        NL.gpu_skin!(w,ex.positions,ex.identities,ex.labels)
        offsets=w.cache.offsets; neighbors=w.cache.neighbors
        stamp=gpu_stage_stamp!(diagnostics,w.backend,:neighbors,stamp)
        if !isempty(ex.centers)
            NL.gpu_inverse_permutation!(w.backend,64)(ex.inverse,w.permutation;ndrange=n)
            NL.gpu_select_indices!(w.backend,64)(ex.sortedcenters,ex.centers,ex.inverse;ndrange=length(ex.centers))
        end
        if !isempty(ex.weights)
            ace_weights!(w.backend,64)(ex.weights,t,ex.parameters;ndrange=length(ex.weights))
        end
        stamp=gpu_stage_stamp!(diagnostics,w.backend,:readout,stamp)
        for firstsite=1:ex.batch_size:length(ex.centers)
            count=min(ex.batch_size,length(ex.centers)-firstsite+1)
            args=(w.x,w.species,w.permutation,w.lengths,w.periodic,w.radius,offsets,neighbors,t,ex.parameters,ex.sortedcenters,firstsite,count)
            # Centers are original indices; convert below before stage dispatch.
            nd=size(ex.U,1)
            if nd>0
                if ex.strategy===:baseline
                    ace_density_baseline!(w.backend,64)(ex.U,ex.V,args...;ndrange=count)
                else
                    ace_density!(w.backend,64)(ex.U,ex.V,args...;ndrange=64*count)
                end
            end
            stamp=gpu_stage_stamp!(diagnostics,w.backend,:density,stamp)
            if ex.strategy===:compact
                ace_compact!(w.backend,64)(w.force,w.partials,ex.dU,ex.dV,ex.U,ex.V,ex.weights,args...,energy,forces,virial,Val(t.degree);ndrange=cld(count,4)*64)
                stamp=gpu_stage_stamp!(diagnostics,w.backend,:fused_products_edges,stamp)
            else
                ace_products!(w.backend,64)(w.partials,ex.dU,ex.dV,ex.U,ex.V,ex.weights,w.species,t,ex.sortedcenters,firstsite,count,forces||virial,energy,Val(t.degree);ndrange=count)
                stamp=gpu_stage_stamp!(diagnostics,w.backend,:products,stamp)
                ace_edges!(w.backend,64)(w.force,w.partials,ex.dU,ex.dV,args...,energy,forces,virial;ndrange=count)
                stamp=gpu_stage_stamp!(diagnostics,w.backend,:edges,stamp)
            end
        end
        NL.gpu_reduce!(w.totals,w.partials,w.reducetemp)
        stamp=gpu_stage_stamp!(diagnostics,w.backend,:reduction,stamp)
        status=only(Array(view(w.totals,1:1,1:1)))
        isfinite(status) || throw(DomainError(status,"nonfinite GPU ACE or coincident atoms"))
    end
        own(a)=owning ? copy(a) : a
    merge(energy ? (;energy=own(view(w.totals,1:1,:))) : (;),forces ? (;forces=own(w.force)) : (;),
          virial ? (;virial=reshape(own(view(w.totals,2:10,:)),3,3)) : (;))
end
gpu_result(y::NamedTuple)=map(Array,y)

@inline gpu_radial(::Val{:mixed},ch,r,n,iz,jz,ps)=tuple_radial(ch.channels,r,n,iz+(jz-1)*ch.nz,ps)
@inline function tuple_radial(channels::Tuple,r,n,index,ps)
    if index==1
        ch=first(channels); gpu_radial(ch.kind,ch,r,n,1,1,ps)
    else
        tuple_radial(Base.tail(channels),r,n,index-1,ps)
    end
end
@inline tuple_radial(::Tuple{},r,n,index,ps)=(typeof(r)(NaN),typeof(r)(NaN))
@inline gpu_radial(::Val{:spline_single},ch,r,n,iz,jz,ps)=gpu_radial(ch.data,r,n,1,1,ps)
@inline gpu_radial(::Val{:subset},ch,r,n,iz,jz,ps)=gpu_radial(ch.parent,r,(@inbounds ch.indices[n]),1,1,ps)
@inline function gpu_callable(ch,r,ps,::Val{N}) where N
    T=typeof(r); v=MVector{N,T}(undef); d=MVector{N,T}(undef)
    ch.evaluator(v,d,r,ps,ch.config)
    v,d
end
@inline function gpu_radial(::Val{:callable},ch,r,n,iz,jz,ps)
    r>=ch.cutoff && return (zero(r),zero(r))
    v,d=gpu_callable(ch,r,ps,ch.count)
    @inbounds (v[n],d[n])
end
@inline function gpu_radial(::Val{:tabulated},ch,r,n,iz,jz,ps)
    r>=ch.cutoff && return (zero(r),zero(r))
    @inbounds begin
        r<ch.nodes[1] && return (typeof(r)(NaN),typeof(r)(NaN))
        lo=1; hi=length(ch.nodes)
        while lo<hi
            mid=(lo+hi+1)>>>1
            if ch.nodes[mid]<=r; lo=mid; else; hi=mid-1; end
        end
        k=min(lo,length(ch.nodes)-1); h=ch.nodes[k+1]-ch.nodes[k]; t=(r-ch.nodes[k])/h
        c=ntuple(j->ch.coefficients[j,n,k],4)
        (((c[4]*t+c[3])*t+c[2])*t+c[1],(3*c[4]*t*t+2*c[3]*t+c[2])/h)
    end
end

check_gpu_data(x,backend,T)=nothing
check_gpu_data(x::AbstractFloat,backend,T)=typeof(x)===T ? nothing : throw(ArgumentError("scalar parameter precision mismatch"))
check_gpu_data(x::Tuple,backend,T)=foreach(v->check_gpu_data(v,backend,T),x)
function check_gpu_data(x::AbstractArray,backend,T)
    KA.get_backend(x)==backend || throw(ArgumentError("model/parameter backend mismatch; adapt explicitly"))
    eltype(x)<:AbstractFloat && eltype(x)!==T && throw(ArgumentError("parameter precision mismatch"))
end
check_gpu_data(x::NamedTuple,backend,T)=foreach(v->check_gpu_data(v,backend,T),values(x))
function update_gpu_data!(dst::AbstractArray,src::AbstractArray)
    size(dst)==size(src) || throw(DimensionMismatch("parameter shape changed"))
    copyto!(dst,src)
end
update_gpu_data!(dst::Tuple,src::Tuple)=length(dst)==length(src) ? map(update_gpu_data!,dst,src) : throw(DimensionMismatch("parameter tuple"))
function update_gpu_data!(dst::NamedTuple,src::NamedTuple)
    keys(dst)==keys(src) || throw(ArgumentError("parameter structure changed; reprepare"))
    NamedTuple{keys(dst)}(map(update_gpu_data!,values(dst),values(src)))
end
update_gpu_data!(dst,src)=typeof(dst)===typeof(src) ? src : throw(ArgumentError("scalar parameter type changed; adapt and reprepare"))

"""Explicit upload/update boundary; cell/species/topology changes invalidate cached
neighbors. Existing output copies remain owning. Reuses buffers at identical size;
capacity changes allocate outside numerical kernels. Model/parameters stay resident.
"""
function prepare_executor!(ex::GPUExecutor,sys;centers=nothing)
    w=ex.workspace; T=ex.model.precision
    check_execution_device(ex.model.tables.reference)
    input=NL.gpu_inputs(sys,NI.search_cutoff(ex.model.source),T,w.config)
    for s in input.species; element_channel(ex.model.source,s); end
    n=length(sys); selected=collect(center_indices(centers,n))
    upload(dst,src)=begin
        buffer=size(dst)==size(src) ? dst : KA.allocate(w.backend,eltype(src),size(src)...)
        copyto!(buffer,src); buffer
    end
    ex.positions=upload(ex.positions,input.positions)
    ex.identities=upload(ex.identities,input.identities)
    ex.centers=upload(ex.centers,selected)
    ex.sortedcenters=upload(ex.sortedcenters,selected)
    ex.inverse=upload(ex.inverse,zeros(Int,n))
    w.lengths=input.lengths; w.periodic=input.periodic; w.counts=input.counts; w.radius=input.radius
    NL.gpu_resize!(w,n); w.cache=nothing
    ex
end

mutable struct GPUACECalculator{M,P,S}
    model::M
    parameters::P
    state::S
    executor::Any
end
AC.energy_unit(::GPUACECalculator)=u"eV"
AC.length_unit(::GPUACECalculator)=u"Å"
AC.promote_force_type(sys::AbstractSystem{3},c::GPUACECalculator)=SVector{3,typeof(zero(c.model.precision)*u"eV/Å")}
AC.zero_energy(sys::AbstractSystem{3},c::GPUACECalculator)=zero(c.model.precision)*u"eV"
AC.zero_virial(sys::AbstractSystem{3},c::GPUACECalculator)=zero(SMatrix{3,3,typeof(AC.zero_energy(sys,c)),9})
AC.get_parameters(c::GPUACECalculator)=deepcopy(c.parameters)
function AC.set_parameters!(c::GPUACECalculator,ps)
    check_gpu_data(ps,KA.get_backend(c.model.tables.reference),c.model.precision)
    check_placement(ps,c.model.tables.reference)
    check_update_shapes(c.parameters,ps)
    c.parameters=update_gpu_data!(c.parameters,ps); c
end
function calculate_gpu_ace(props,sys,c;parameters=c.parameters,state=c.state,centers=nothing,kwargs...)
    props=NL.normalize_properties(props isa Tuple ? props : (props,))
    isempty(props) && return (;state)
    if isnothing(c.executor) || !isempty(kwargs)
        c.executor=GPUExecutor(sys,c.model,parameters;centers,kwargs...)
    else
        prepare_executor!(c.executor,sys;centers)
    end
    r=gpu_result(execute_gpu!(c.executor;parameters,properties=props))
    merge(haskey(r,:energy) ? (;energy=only(r.energy)*u"eV") : (;),
          haskey(r,:forces) ? (;forces=force_vectors(r.forces*u"eV/Å")) : (;),
          haskey(r,:virial) ? (;virial=static_virial(r.virial*u"eV")) : (;),(;state))
end
AC.calculate(props::Union{Tuple,AC.Energy,AC.Forces,AC.Virial},sys::AbstractSystem{3},c::GPUACECalculator;kwargs...)=calculate_gpu_ace(props,sys,c;kwargs...)
AC.potential_energy(sys::AbstractSystem{3},c::GPUACECalculator;kwargs...)=calculate_gpu_ace((AC.Energy(),),sys,c;kwargs...).energy
AC.forces(sys::AbstractSystem{3},c::GPUACECalculator;kwargs...)=calculate_gpu_ace((AC.Forces(),),sys,c;kwargs...).forces
AC.virial(sys::AbstractSystem{3},c::GPUACECalculator;kwargs...)=calculate_gpu_ace((AC.Virial(),),sys,c;kwargs...).virial
AC.energy_forces_virial(sys::AbstractSystem{3},c::GPUACECalculator;kwargs...)=Base.structdiff(calculate_gpu_ace((AC.Energy(),AC.Forces(),AC.Virial()),sys,c;kwargs...),NamedTuple{(:state,)})

# Experimental partial fusion: four 16-lane site groups per workgroup. An owner
# contracts each site's products, then a workgroup barrier publishes global
# adjoints. All sixteen site lanes revisit edges and cooperatively reduce central
# force/virial/pair energy. Only 13*4*16 scalars occupy shared memory. No hardware
# warp/subgroup width is presumed; KA barriers cover the complete workgroup.
@kernel function ace_compact!(force,partials,dU,dV,U,V,weights,x,channels,perm,L,pb,c,offsets,neighbors,tables,ps,centers,firstsite,count,energy,forces,virial,::Val{D}) where D
    state=@private Int (3,) # site in group, site lane, batch slot
    group=@index(Group,Linear); lane=@index(Local,Linear)
    state[1]=(lane-1)÷16+1; state[2]=(lane-1)%16+1; state[3]=(group-1)*4+state[1]
    sums=@localmem eltype(x) (13,4,16)
    @inbounds begin
        if state[3]<=count && state[2]==1
            gpu_site_products!(state[3],partials,dU,dV,U,V,weights,channels,tables,centers,firstsite,count,forces||virial,energy,Val(D))
        end
        @synchronize
        values=ntuple(_->zero(eltype(x)),Val(13))
        if state[3]<=count
            i=centers[firstsite+state[3]-1]; iz=channels[i]
            firstedge=(i==1 ? 0 : offsets[i-1])+state[2]
            values=gpu_edge_partial(force,dU,dV,x,channels,perm,L,pb,c,offsets,neighbors,tables,ps,i,iz,state[3],firstedge,16,energy,forces,virial)
        end
        for a=1:13; sums[a,state[1],state[2]]=values[a]; end
        @synchronize
        for step in (8,4,2,1)
            if state[2]<=step
                for a=1:13; sums[a,state[1],state[2]]+=sums[a,state[1],state[2]+step]; end
            end
            @synchronize
        end
        if state[3]<=count && state[2]==1
            i=centers[firstsite+state[3]-1]
            partials[1,i]+=sums[1,state[1],1]
            if forces; for a=1:3; @atomic force[a,perm[i]]+=sums[a+1,state[1],1]; end; end
            if virial; for a=1:9; partials[a+1,i]=sums[a+4,state[1],1]; end; end
        end
    end
end

# Array payload accounting, distinct from host allocation measurements.
array_payload(x::AbstractArray)=sizeof(x)
array_payload(x::NamedTuple)=sum(array_payload,values(x);init=0)
array_payload(x::Tuple)=sum(array_payload,x;init=0)
array_payload(x)=0
function gpu_memory(ex::GPUExecutor)
    w=ex.workspace; cache=w.cache
    (;model_tables=array_payload(ex.model.tables),parameters=array_payload(ex.parameters),
      host_model_descriptor=Base.summarysize(ex.model.source),inputs=sizeof(ex.positions)+sizeof(ex.identities)+sizeof(ex.labels),
      density=sizeof(ex.U)+sizeof(ex.V),adjoints=sizeof(ex.dU)+sizeof(ex.dV),weights=sizeof(ex.weights),
      outputs=sizeof(w.force)+sizeof(w.partials)+sizeof(w.totals),
      connectivity=isnothing(cache) ? 0 : sizeof(cache.offsets)+sizeof(cache.neighbors),
      sorted_cells=sizeof(w.records)+sizeof(w.x)+sizeof(w.ends)+sizeof(w.species)+sizeof(w.permutation),
      selection=sizeof(ex.centers)+sizeof(ex.sortedcenters)+sizeof(ex.inverse),
      note="Array payloads; sorting/reduction/skin auxiliary allocations and host object headers are additional")
end

@inline gpu_stage_stamp!(::Nothing,backend,name,stamp)=stamp
function gpu_stage_stamp!(diagnostics::AbstractDict,backend,name,stamp)
    KA.synchronize(backend); current=time_ns()
    diagnostics[name]=get(diagnostics,name,0.)+(current-stamp)*1e-9
    current
end

# Optional vendor extensions refine generic backend identity with the actual device.
storage_device(a)=KA.get_backend(a)
check_execution_device(a)=nothing
check_placement(x,reference)=nothing
check_placement(x::NamedTuple,reference)=foreach(v->check_placement(v,reference),values(x))
check_placement(x::Tuple,reference)=foreach(v->check_placement(v,reference),x)
function check_placement(x::AbstractArray,reference)
    storage_device(x)==storage_device(reference) || throw(ArgumentError("model, parameters and scratch must be on the same selected device; adapt explicitly"))
end

AC.calculate(props::Union{Tuple,AC.Energy,AC.Forces,AC.Virial},sys::AbstractSystem{3},c::GPUACECalculator,parameters;state=c.state,kwargs...)=calculate_gpu_ace(props,sys,c;parameters,state,kwargs...)
AC.calculate(props::Union{Tuple,AC.Energy,AC.Forces,AC.Virial},sys::AbstractSystem{3},c::GPUACECalculator,parameters,state;kwargs...)=calculate_gpu_ace(props,sys,c;parameters,state,kwargs...)
function AC.energy_forces(sys::AbstractSystem{3},c::GPUACECalculator;kwargs...)
    r=calculate_gpu_ace((AC.Energy(),AC.Forces()),sys,c;kwargs...); (;energy=r.energy,forces=r.forces)
end
function AC.forces!(f::AbstractVector,sys::AbstractSystem{3},c::GPUACECalculator;kwargs...)
    check_destination(f,sys); add_forces!(f,AC.forces(sys,c;kwargs...))
end
function AC.energy_forces!(f::AbstractVector,sys::AbstractSystem{3},c::GPUACECalculator;kwargs...)
    check_destination(f,sys); r=AC.energy_forces(sys,c;kwargs...)
    add_forces!(f,r.forces); (;energy=r.energy,forces=f)
end
function AC.energy_forces_virial!(f::AbstractVector,sys::AbstractSystem{3},c::GPUACECalculator;kwargs...)
    check_destination(f,sys); r=AC.energy_forces_virial(sys,c;kwargs...)
    add_forces!(f,r.forces); (;energy=r.energy,forces=f,virial=r.virial)
end

check_update_shapes(dst::AbstractArray,src::AbstractArray)=size(dst)==size(src) ? nothing : throw(DimensionMismatch("parameter shape changed"))
function check_update_shapes(dst::NamedTuple,src::NamedTuple)
    keys(dst)==keys(src) || throw(ArgumentError("parameter structure changed"))
    foreach(check_update_shapes,values(dst),values(src))
end
function check_update_shapes(dst::Tuple,src::Tuple)
    length(dst)==length(src) || throw(DimensionMismatch("parameter tuple changed"))
    foreach(check_update_shapes,dst,src)
end
check_update_shapes(dst,src)=typeof(dst)===typeof(src) ? nothing : throw(ArgumentError("scalar parameter type changed; adapt and reprepare"))
function save_calculator(path::AbstractString,c::GPUACECalculator;parameters=c.parameters,state=c.state)
    b=adapt_model(c.model,parameters;array_type=Array,precision=c.model.precision)
    cpu=ACECalculator(b.model,b.parameters,host_data(state);precision=c.model.precision,angular_backend=:block,slots=1)
    save_calculator(path,cpu)
end
storage_device(a::SubArray)=storage_device(parent(a))
storage_device(a::Base.ReshapedArray)=storage_device(parent(a))
storage_device(a::Base.ReinterpretArray)=storage_device(parent(a))
