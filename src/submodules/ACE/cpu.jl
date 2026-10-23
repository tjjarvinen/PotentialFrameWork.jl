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
end
function SiteScratch(model,capacity,::Type{T}) where T
    nr=length(model.radial_spec); nq=length(model.pair_spec); ny=(model.lmax+1)^2
    SiteScratch(zeros(SVector{3,T},capacity),zeros(T,capacity,nr),zeros(T,capacity,nr),
        zeros(T,capacity,nq),zeros(T,capacity,nq),zeros(T,capacity,ny),zeros(SVector{3,T},capacity,ny),
        zeros(T,length(model.density)),zeros(T,length(model.products)),zeros(T,length(model.density)))
end
mutable struct CPUExecutor{T,B,M}
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
end
function CPUExecutor(model,env;precision=Float64,slots=Threads.nthreads(),angular_backend=:sphericart)
    precision in (Float32,Float64) || throw(ArgumentError("CPU executor precision must be Float32/64"))
    slots>0 || throw(ArgumentError("slots must be positive"))
    basis=make_angular(model,precision,angular_backend,env.width)
    T=precision
    w=CPUExecutor(model,env,basis,[SiteScratch(model,0,T) for _=1:slots],zeros(T,3,0),Int[],Int[],
        zeros(T,length(model.products),length(model.atomic_numbers)),T[],zeros(T,9,0),0,UInt(0))
    prepare_executor!(w,env)
end
function prepare_executor!(w::CPUExecutor{T},env) where T
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
    size(w.virials,2)==n || (w.virials=zeros(T,9,n))
    w.environment_version=env.version
    w
end

function angular_ed! end

# SIMD Horner evaluation of the actual resolved cubic spline. Transform and
# interval selection are scalar per lane; channel/species selection is outside
# SIMD arithmetic. This is intentionally separate from an analytic recurrence
# for a polynomial basis that is not the reconstructed model.
function radial_block!(R,dR,rs,ch,start,::Val{W},active) where W
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
function project_rows!(s,model,lo,hi)
    for (a,(r,m)) in enumerate(model.density),j=lo:hi
        s.A[a]+=s.R[j,r]*s.Y[j,m]
    end
end
function cpu_site!(w::CPUExecutor{T},s::SiteScratch,i,X,H,ps,::Val{W},simd,fused=false) where {T,W}
    env=w.environment; model=w.model; iz=env.species[i]
    firstedge=env.offsets[1,i]; count=env.offsets[end,i]-firstedge
    for j=1:count
        k=firstedge+j-1; s.Rs[j]=edge_displacement(X,H,env,i,k)
        env.active[k] && iszero(norm(s.Rs[j])) && throw(DomainError(0,"coincident atoms"))
    end
    fill!(s.A,0); fill!(s.dA,0)
    for z in eachindex(model.atomic_numbers)
        lo=env.offsets[z,i]; hi=env.offsets[z+1,i]-1
        if simd
            for k=lo:W:hi
                localstart=k-firstedge+1
                rs=ntuple(lane -> norm(s.Rs[localstart+lane-1]),W)
                active=ntuple(lane -> env.active[k+lane-1],W)
                radial_block!(s.R,s.dR,rs,model.radial[iz,z],localstart,Val(W),active)
                radial_block!(s.Q,s.dQ,rs,model.pair_radial[iz,z],localstart,Val(W),active)
                if fused
                    stop=localstart+W-1
                    scale_rows!(s,ps,iz,z,localstart,stop)
                    angular_ed!(view(s.Y,localstart:stop,:),view(s.dY,localstart:stop,:),w.basis,view(s.Rs,localstart:stop))
                    project_rows!(s,model,localstart,stop)
                end
            end
        else
            for k=lo:hi
                j=k-firstedge+1
                if env.active[k]
                    s.R[j,:],s.dR[j,:]=radial_ed(model.radial[iz,z],norm(s.Rs[j]))
                    s.Q[j,:],s.dQ[j,:]=radial_ed(model.pair_radial[iz,z],norm(s.Rs[j]))
                else
                    for a=1:size(s.R,2); s.R[j,a]=s.dR[j,a]=0; end
                    for a=1:size(s.Q,2); s.Q[j,a]=s.dQ[j,a]=0; end
                end
            end
        end
        fused || scale_rows!(s,ps,iz,z,lo-firstedge+1,hi-firstedge+1)
    end
    if !fused && count>0
        angular_ed!(view(s.Y,1:count,:),view(s.dY,1:count,:),w.basis,view(s.Rs,1:count))
    end
    fused || project_rows!(s,model,1,count)
    energy=T(model.reference[iz])
    for (k,indices) in enumerate(model.products)
        value=one(T)
        for a in indices; value*=s.A[a]; end
        s.P[k]=value; energy+=w.weights[k,iz]*value
        for t in eachindex(indices)
            v=w.weights[k,iz]
            for u in eachindex(indices); u==t || (v*=s.A[indices[u]]); end
            s.dA[indices[t]]+=v
        end
    end
    for j=1:count,n=1:size(s.Q,2); energy+=s.Q[j,n]*ps.Wpair[n,iz]; end
    w.energies[i]=energy
    for j=1:count
        k=firstedge+j-1; env.active[k] || continue
        d=s.Rs[j]; direction=d/norm(d); g=zero(SVector{3,T})
        for (a,(r,m)) in enumerate(model.density)
            g+=s.dA[a]*(s.dR[j,r]*s.Y[j,m]*direction+s.R[j,r]*s.dY[j,m])
        end
        for n=1:size(s.Q,2); g+=ps.Wpair[n,iz]*s.dQ[j,n]*direction; end
        dest=env.neighbours[k]
        for a=1:3
            w.force[a,w.inverse[i]]+=g[a]; w.force[a,w.inverse[dest]]-=g[a]
            for b=1:3; w.virials[a+3(b-1),i]-=d[a]*g[b]; end
        end
    end
    nothing
end

"""Prepared CPU inference with original-order numeric inputs and one sorted
shared force array. `prepare_executor!` follows every environment rebuild.
Precontraction is refreshed from authoritative parameters on every execution.
Scratch is exclusive per slot; phases advertise complete endpoint footprints.
"""
function execute_prepared!(w::CPUExecutor,X,H,ps;parallel=false,simd=true,fused=false)
    fused && !simd && throw(ArgumentError("tiled fusion requires SIMD radial blocks"))
    env=w.environment
    env.version==w.environment_version || throw(ArgumentError("prepare executor after rebuilding the environment"))
    length(env.species)==length(w.p) || throw(ArgumentError("prepare executor after environment changes"))
    maximum((env.offsets[end,i]-env.offsets[1,i] for i in eachindex(env.species));init=0)<=w.capacity ||
        throw(ArgumentError("prepare executor after environment growth"))
    mul!(w.weights,transpose(w.model.coupling),ps.WB)
    fill!(w.force,0); fill!(w.energies,0); fill!(w.virials,0)
    if parallel
        for phase in env.phases
            slots=min(length(w.scratch),length(phase))
            if slots<=1
                for i in phase
                    cpu_site!(w,w.scratch[1],i,X,H,ps,Val(env.width),simd,fused)
                end
            else
                Threads.@threads :dynamic for slot=1:slots
                    for k=slot:slots:length(phase)
                        cpu_site!(w,w.scratch[slot],phase[k],X,H,ps,Val(env.width),simd,fused)
                    end
                end
            end
        end
    else
        for i in eachindex(env.species); cpu_site!(w,w.scratch[1],i,X,H,ps,Val(env.width),simd,fused); end
    end
    w
end
"""Owning original-order numerical outputs. Prepared scratch is borrowed only
inside the executor; this result survives subsequent executions and resizing.
"""
function prepared_result(w::CPUExecutor)
    (;energy=sum(w.energies),forces=w.force[:,w.inverse],virial=reshape(vec(sum(w.virials;dims=2)),3,3))
end
