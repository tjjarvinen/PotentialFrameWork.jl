# Scratch belongs to an exclusive task slot. Neither ownership nor indexing uses
# Threads.threadid(); a task may migrate while retaining its slot.
mutable struct SiteScratch{T}
    Rs::Vector{SVector{3,T}}
    R::Matrix{T}
    dR::Matrix{T}
    Q::Matrix{T}
    dQ::Matrix{T}
    Y::Matrix{T}
    dY::Matrix{SVector{3,T}}
    A::Vector{T}
    P::Vector{T}
    dA::Vector{T}
    U::Vector{T}
    V::Vector{T}
    dV::Vector{T}
    polynomial::Matrix{T}
    reverse_polynomial::Matrix{T}
end
function SiteScratch(model,capacity,::Type{T}) where T
    nr=length(model.radial_spec); nq=length(model.pair_spec); ny=(model.lmax+1)^2
    SiteScratch(zeros(SVector{3,T},capacity),zeros(T,capacity,nr),zeros(T,capacity,nr),
        zeros(T,capacity,nq),zeros(T,capacity,nq),zeros(T,capacity,ny),zeros(SVector{3,T},capacity,ny),
        zeros(T,length(model.density)),zeros(T,length(model.products)),zeros(T,length(model.density)),zeros(T,length(model.density)),zeros(T,length(model.density)),
        zeros(T,length(model.density)),projection_work(model,T)...)
end
struct PreparedACE{M,B,C,W}
    model::M
    basis::B
    weights::C
    simd::Bool
    fused::Bool
end
mutable struct CPUExecutor{T,B,M,R,S,P}
    model::M
    environment::FullEnvironment
    basis::B
    scratch::Vector{SiteScratch{T}}
    force::Matrix{T} # the only shared numerical force array, sorted by p
    p::Vector{Int}
    inverse::Vector{Int}
    weights::Matrix{T}
    energies::Vector{T}
    virials::Matrix{T}
    capacity::Int
    environment_version::UInt
    readout::R
    center_mask::BitVector
    slots::S
    batch_size::Int
    prepared::P
end
function CPUExecutor(model,env;precision=Float64,slots=Threads.nthreads(),angular_backend=:sphericart,batch_size=1)
    precision in (Float32,Float64) || throw(ArgumentError("CPU executor precision must be Float32/64"))
    slots>0 || throw(ArgumentError("slots must be positive"))
    batch_size isa Integer && batch_size>0 || throw(ArgumentError("batch_size must be positive"))
    NL.check_units(model)
    basis=make_angular(model,precision,angular_backend,env.width)
    T=precision
    scratch=[NI.make_scratch(model,0,T) for _=1:slots]
    weights=zeros(T,length(model.products),length(model.atomic_numbers))
    prepared=PreparedACE{typeof(model),typeof(basis),typeof(weights),env.width}(model,basis,weights,true,false)
    w=CPUExecutor(model,env,basis,scratch,zeros(T,3,0),Int[],Int[],
        weights,T[],zeros(T,9,0),0,UInt(0),
        zeros(promote_type(Float64,eltype(model.coupling)),size(model.coupling,1),length(model.atomic_numbers)),BitVector(),NL.site_slots(scratch,0,batch_size,T),Int(batch_size),prepared)
    prepare_executor!(w,env)
end
function prepare_executor!(w::CPUExecutor{T},env) where T
    env.width==prepared_width(w.prepared) || throw(ArgumentError("prepare a new executor for a changed layout width"))
    w.environment=env; n=length(env.species)
    capacity=maximum((env.offsets[end,i]-env.offsets[1,i] for i=1:n);init=0)
    if capacity>w.capacity
        w.capacity=max(capacity,2w.capacity)
        for k in eachindex(w.scratch); w.scratch[k]=SiteScratch(w.model,w.capacity,T); end
    end
    w.p=sortperm(env.species); resize!(w.inverse,n)
    for k in eachindex(w.p); w.inverse[w.p[k]]=k; end
    size(w.force,2)==n || (w.force=zeros(T,3,n))
    resize!(w.energies,n)
    resize!(w.center_mask,n)
    size(w.virials,2)==n || (w.virials=zeros(T,9,n))
    for k in eachindex(w.slots)
        w.slots[k].scratch=w.scratch[k]
        for out in w.slots[k].outputs; resize!(out.gradients,w.capacity); end
    end
    w.environment_version=env.version
    w
end

function angular_ed! end

# SIMD Horner evaluation of the actual resolved cubic spline. Transform and
# interval selection are scalar per lane; channel/species selection is outside
# SIMD arithmetic. This is intentionally separate from an analytic recurrence
# for a polynomial basis that is not the reconstructed model.
function radial_block!(R,dR,rs,ch::RadialChannel,start,::Val{W},active,ps=NamedTuple()) where W
    xs=ntuple(k -> transformed(ch.transform,rs[k]),W)
    idx=ntuple(k -> clamp(searchsortedlast(ch.nodes,xs[k]),1,length(ch.nodes)-1),W)
    us=ntuple(k -> (xs[k]-ch.nodes[idx[k]])/(ch.nodes[idx[k]+1]-ch.nodes[idx[k]]),W)
    # Scalar joint transform/envelope derivative samples. Polynomial arithmetic
    # below is SIMD; scalar derivative path remains the transparent baseline.
    factors=ntuple(W) do k
        t=ch.transform; r=rs[k]; x=xs[k]; s=(r-t.rin)/(t.r0-t.rin); h=t.q-t.p
        dx=if r<=t.rin || x<=-1 || x>=1
            zero(r)
        else
            den=1+s^h; transformed_u=s^t.q/den
            transform_deriv=(t.q*s^(t.q-1)/den-(h==0 ? zero(s) : h*s^(t.q+h-1)/den^2))/(t.r0-t.rin)
            -2t.a*transform_deriv/((1+t.a*transformed_u)^2*(t.ycut-t.yin))
        end
        envelope_val=envelope(ch,r,x); ep=ch.envelope_parameters
        envelope_deriv=if ch.envelope===:PolyEnvelope2sX
            p=Int(ep[3]); q=Int(ep[4])
            ep[1]<x<ep[2] ? ep[5]*(p*(x-ep[1])^(p-1)*(ep[2]-x)^q-q*(x-ep[1])^p*(ep[2]-x)^(q-1))*dx : zero(x)
        else
            rc,r0,p=ep; p=Int(p)
            r>rc ? zero(r) : p/r0*((rc/r0)^(-p-1)-(r/r0)^(-p-1))
        end
        (envelope_val,envelope_deriv,dx/(ch.nodes[idx[k]+1]-ch.nodes[idx[k]]))
    end
    u=Vec(us); e=Vec(ntuple(k -> factors[k][1],W)); de=Vec(ntuple(k -> factors[k][2],W))
    du=Vec(ntuple(k -> factors[k][3],W))
    for n=1:size(R,2)
        c=ntuple(a -> Vec(ntuple(k -> ch.coefficients[a,n,idx[k]],W)),4)
        v=((c[4]*u+c[3])*u+c[2])*u+c[1]
        dv=(3c[4]*u*u+2c[3]*u+c[2])*du
        value=v*e; derivative=dv*e+v*de
        for k=1:W
            R[start+k-1,n]=active[k] ? value[k] : zero(eltype(R))
            dR[start+k-1,n]=active[k] ? derivative[k] : zero(eltype(dR))
        end
    end
end

function scale_rows!(s,ps,iz,z,lo,hi)
    for j=lo:hi
        for n=1:size(s.R,2)
            scale=radial_scale(ps,:radial,n,iz,z)
            s.R[j,n]*=scale; s.dR[j,n]*=scale
        end
        for n=1:size(s.Q,2)
            scale=radial_scale(ps,:pair,n,iz,z)
            s.Q[j,n]*=scale; s.dQ[j,n]*=scale
        end
    end
end
function project_rows!(s,model,lo,hi,selected=false)
    for (a,(r,m)) in enumerate(model.density),j=lo:hi
        value=s.R[j,r]*s.Y[j,m]; s.A[a]+=value
        if selected; s.V[a]+=value; else; s.U[a]+=value; end
    end
end
function NI.evaluate_site!(out,p::PreparedACE{M,B,C,W},site,ps,s::SiteScratch{T},props) where {M,B,C,W,T}
    model=p.model; iz=NI.central_species(site)
    count=NI.neighbour_count(site); firstedge=1
    simd=p.simd; fused=p.fused
    for j=1:count
        k=firstedge+j-1; s.Rs[j]=NI.displacement(site,j)
        NI.active(site,k) && iszero(norm(s.Rs[j])) && throw(DomainError(0,"coincident atoms"))
    end
    fill!(s.A,0); fill!(s.U,0); fill!(s.V,0); fill!(s.dA,0)
    for z in eachindex(model.atomic_numbers)
        range=NI.species_range(site,z); lo=first(range); hi=last(range)
        if simd
            for k=lo:W:hi
                localstart=k-firstedge+1
                rs=ntuple(lane -> norm(s.Rs[localstart+lane-1]),W)
                active=ntuple(lane -> NI.active(site,k+lane-1),W)
                radial_block!(s.R,s.dR,rs,model.radial[iz,z],localstart,Val(W),active,radial_parameters(ps,:custom_radial))
                radial_block!(s.Q,s.dQ,rs,model.pair_radial[iz,z],localstart,Val(W),active,radial_parameters(ps,:custom_pair))
                if fused
                    stop=localstart+W-1
                    scale_rows!(s,ps,iz,z,localstart,stop)
                    angular_ed!(view(s.Y,localstart:stop,:),view(s.dY,localstart:stop,:),p.basis,view(s.Rs,localstart:stop))
                    project_rows!(s,model,localstart,stop,z in model.policy.species)
                end
            end
        else
            for k=lo:hi
                j=k-firstedge+1
                if NI.active(site,k)
                    s.R[j,:],s.dR[j,:]=radial_ed(model.radial[iz,z],norm(s.Rs[j]),radial_parameters(ps,:custom_radial))
                    s.Q[j,:],s.dQ[j,:]=radial_ed(model.pair_radial[iz,z],norm(s.Rs[j]),radial_parameters(ps,:custom_pair))
                else
                    for a=1:size(s.R,2); s.R[j,a]=s.dR[j,a]=0; end
                    for a=1:size(s.Q,2); s.Q[j,a]=s.dQ[j,a]=0; end
                end
            end
        end
        fused || scale_rows!(s,ps,iz,z,lo-firstedge+1,hi-firstedge+1)
        if !pair_allowed(model,z)
            for j=lo-firstedge+1:hi-firstedge+1,n=1:size(s.Q,2)
                s.Q[j,n]=s.dQ[j,n]=0
            end
        end
    end
    if !fused && count>0
        angular_ed!(view(s.Y,1:count,:),view(s.dY,1:count,:),p.basis,view(s.Rs,1:count))
    end
    if !fused
        for z in eachindex(model.atomic_numbers)
            project_rows!(s,model,first(NI.species_range(site,z)),last(NI.species_range(site,z)),z in model.policy.species)
        end
    end
    energy=T(reference_value(model,iz))
    product_projection!(s.P,s.dA,s.dV,s.U,s.V,model.products,view(p.weights,:,iz),model.policy.orders,s.polynomial,s.reverse_polynomial)
    for k in eachindex(s.P); energy+=p.weights[k,iz]*s.P[k]; end
    for j=1:count,n=1:size(s.Q,2); energy+=s.Q[j,n]*ps.Wpair[n,iz]; end
    NL.wants(AC.Energy,props) && NI.set_energy!(out,energy)
    (NL.wants(AC.Forces,props) || NL.wants(AC.Virial,props)) || return nothing
    for j=1:count
        k=firstedge+j-1; NI.active(site,k) || continue
        d=s.Rs[j]; direction=d/norm(d); g=zero(SVector{3,T})
        for (a,(r,m)) in enumerate(model.density)
            adj=NI.neighbour_species(site,k) in model.policy.species ? s.dV : s.dA
            g+=adj[a]*(s.dR[j,r]*s.Y[j,m]*direction+s.R[j,r]*s.dY[j,m])
        end
        for n=1:size(s.Q,2); g+=ps.Wpair[n,iz]*s.dQ[j,n]*direction; end
        NI.set_gradient!(out,j,g)
    end
    nothing
end

# Compatibility helper for existing site-footprint diagnostics.
function cpu_site!(w::CPUExecutor,s::SiteScratch,i,X,H,ps,::Val{W},simd,fused=false) where W
    p=PreparedACE{typeof(w.model),typeof(w.basis),typeof(w.weights),W}(w.model,w.basis,w.weights,simd,fused)
    site=NL.SiteView(w.environment,X,H,i); out=w.slots[1].outputs[1]
    props=(AC.Energy(),AC.Forces(),AC.Virial())
    NI.evaluate_site!(out,p,site,ps,s,props)
    NL.reduce_site!(w.force,w.inverse,w.energies,w.virials,w.environment,i,site,out,props)
end

"""Prepared CPU inference with original-order numeric inputs and one sorted
shared force array. `prepare_executor!` follows every environment rebuild.
Precontraction is refreshed from authoritative parameters on every execution.
Scratch is exclusive per slot; phases advertise complete endpoint footprints.
"""
function execute_prepared_impl!(w::CPUExecutor,X,H,ps;parallel=false,simd=true,fused=false,centers=nothing,properties=(AC.Energy(),AC.Forces(),AC.Virial()))
    fused && !simd && throw(ArgumentError("tiled fusion requires SIMD radial blocks"))
    env=w.environment
    env.version==w.environment_version || throw(ArgumentError("prepare executor after rebuilding the environment"))
    length(env.species)==length(w.p) || throw(ArgumentError("prepare executor after environment changes"))
    maximum((env.offsets[end,i]-env.offsets[1,i] for i in eachindex(env.species));init=0)<=w.capacity ||
        throw(ArgumentError("prepare executor after environment growth"))
    selected=center_indices(centers,length(env.species))
    selected_mask=w.center_mask; fill!(selected_mask,false)
    for i in selected; selected_mask[i]=true; end
    size(ps.WB)==size(w.readout) || throw(DimensionMismatch("CPU readout shape"))
    for z in axes(w.readout,2),i in axes(w.readout,1)
        w.readout[i,z]=w.model.policy.mask[i] ? ps.WB[i,z] : zero(eltype(w.readout))
    end
    mul!(w.weights,transpose(w.model.coupling),w.readout)
    props=NL.normalize_properties(properties)
    env.width==prepared_width(w.prepared) || throw(ArgumentError("stale layout width"))
    p=typeof(w.prepared)(w.model,w.basis,w.weights,simd,fused)
    NL.execute_sites!(w.force,w.inverse,w.energies,w.virials,env,selected,selected_mask,w.slots,p,X,H,ps,props,w.batch_size;parallel)

    w
end
"""Owning original-order numerical outputs. Prepared scratch is borrowed only
inside the executor; this result survives subsequent executions and resizing.
"""
function prepared_result(w::CPUExecutor)
    (;energy=sum(w.energies),forces=w.force[:,w.inverse],virial=reshape(vec(sum(w.virials;dims=2)),3,3))
end

# Callable adapter: scalar lanes are explicit, with no claim of arbitrary SIMD AD.
function radial_block!(R,dR,rs,ch::CallableRadial,start,::Val{W},active,ps=NamedTuple()) where W
    for k=1:W
        row=start+k-1
        if active[k]
            radial_evaluate!(view(R,row,:),view(dR,row,:),ch,rs[k],ps)
        else
            for n=1:size(R,2); R[row,n]=dR[row,n]=0; end
        end
    end
end
function radial_block!(R,dR,rs,ch,start,::Val{W},active,ps=NamedTuple()) where W
    for k=1:W
        row=start+k-1
        if active[k]
            v,d=radial_ed(ch,rs[k],ps)
            for n in eachindex(v); R[row,n]=v[n]; dR[row,n]=d[n]; end
        else
            for n=1:size(R,2); R[row,n]=dR[row,n]=0; end
        end
    end
end

function execute_prepared!(w::CPUExecutor,X,H,ps;kwargs...)
    has_source(w.model) ? Base.invokelatest(execute_prepared_impl!,w,X,H,ps;kwargs...) : execute_prepared_impl!(w,X,H,ps;kwargs...)
end

AC.energy_unit(::ResolvedModel)=u"eV"
AC.length_unit(::ResolvedModel)=u"Å"
AC.energy_unit(p::PreparedACE)=AC.energy_unit(p.model)
AC.length_unit(p::PreparedACE)=AC.length_unit(p.model)
NI.evaluation_kind(::ResolvedModel)=:site
NI.species_channels(m::ResolvedModel)=ChemicalSpecies.(m.atomic_numbers)
NI.search_cutoff(m::ResolvedModel)=maximum(radial_cutoff(ch) for family in (m.radial,m.pair_radial) for ch in family)
NI.make_scratch(m::ResolvedModel,capacity,::Type{T}) where T=SiteScratch(m,capacity,T)
NI.make_scratch(p::PreparedACE,capacity,::Type{T}) where T=SiteScratch(p.model,capacity,T)
function NI.prepare_potential(m::ResolvedModel,species;precision=Float64,width=1,angular_backend=:sphericart,simd=true,fused=false,kwargs...)
    fused && !simd && throw(ArgumentError("tiled fusion requires SIMD radial blocks"))
    width in (1,2,4,8,16) || throw(ArgumentError("unsupported environment width"))
    NL.check_units(m)
    # Provider species identities are mapped explicitly by prepare_environment!;
    # numerical channel IDs may only be used with that declared ordering.
    all(s -> element_channel(m,s)>0,species) || throw(ArgumentError("species channels"))
    basis=make_angular(m,precision,angular_backend,width)
    weights=zeros(precision,length(m.products),length(m.atomic_numbers))
    PreparedACE{typeof(m),typeof(basis),typeof(weights),width}(m,basis,weights,simd,fused)
end
function NI.update_prepared!(p::PreparedACE,ps)
    readout=ps.WB .* collect(Bool,p.model.policy.mask)
    mul!(p.weights,transpose(p.model.coupling),readout)
    p
end

prepared_width(::PreparedACE{M,B,C,W}) where {M,B,C,W}=W
