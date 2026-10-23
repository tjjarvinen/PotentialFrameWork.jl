function angular_basis(args...)
    throw(ArgumentError("load SpheriCart to use the reference CPU angular backend"))
end
function angular_values end
function angular_ed end

"""Expose nonlinear diagonal radial mixing while preserving the imported readout
coordinates. The resolved spline tables remain fixed; cutoffs are never trainable.
"""
function trainable_parameters(model,ps)
    nz=length(model.atomic_numbers)
    merge(deepcopy(ps),(;radial=(;scale=ones(length(model.radial_spec),nz,nz)),
                       pair=(;scale=ones(length(model.pair_spec),nz,nz))))
end
radial_scale(ps,key,n,iz,jz)=hasproperty(ps,key) ? getproperty(ps,key).scale[n,iz,jz] : 1

function radial_ed(ch::RadialChannel,r::Real)
    t=ch.transform; x=transformed(t,r)
    s=(r-t.rin)/(t.r0-t.rin); h=t.q-t.p
    dx=if r<=t.rin || x<=-1 || x>=1
        zero(r)
    else
        den=1+s^h
        u=s^t.q/den
        du=(t.q*s^(t.q-1)/den-(h==0 ? zero(s) : h*s^(t.q+h-1)/den^2))/(t.r0-t.rin)
        -2t.a*du/((1+t.a*u)^2*(t.ycut-t.yin))
    end
    e=envelope(ch,r,x); ep=ch.envelope_parameters
    de=if ch.envelope===:PolyEnvelope2sX
        if ep[1]<x<ep[2]
            p=Int(ep[3]); q=Int(ep[4])
            ep[5]*(p*(x-ep[1])^(p-1)*(ep[2]-x)^q-q*(x-ep[1])^p*(ep[2]-x)^(q-1))*dx
        else
            zero(x)
        end
    else
        rc,r0,p=ep; p=Int(p)
        r>rc ? zero(r) : p/r0*((rc/r0)^(-p-1)-(r/r0)^(-p-1))
    end
    k=clamp(searchsortedlast(ch.nodes,x),1,length(ch.nodes)-1)
    step=ch.nodes[k+1]-ch.nodes[k]; u=(x-ch.nodes[k])/step; c=ch.coefficients
    values=Vector{typeof(e+u)}(undef,size(c,2)); derivatives=similar(values)
    for n in eachindex(values)
        v=((c[4,n,k]*u+c[3,n,k])*u+c[2,n,k])*u+c[1,n,k]
        dv=(3c[4,n,k]*u^2+2c[3,n,k]*u+c[2,n,k])/step
        values[n]=e*v; derivatives[n]=de*v+e*dv*dx
    end
    values,derivatives
end

function site_diagnostics(X,model,ps,env,H,i;basis=angular_basis(model.lmax,eltype(model.reference)))
    iz=env.species[i]
    edges=[k for k=env.offsets[1,i]:env.offsets[end,i]-1 if env.active[k]]
    Rs=[edge_displacement(X,H,env,i,k) for k in edges]
    zs=[env.species[env.neighbours[k]] for k in edges]
    merge(site_data(model,ps,iz,zs,Rs;basis),(;edges))
end

function site_data(model,ps,iz,zs,Rs;basis=angular_basis(model.lmax,eltype(model.reference)))
    rs=norm.(Rs)
    all(>(0),rs) || throw(DomainError(0,"coincident physical atoms have undefined imported harmonics/pair envelope"))
    T=typeof(zero(eltype(eltype(Rs)))+parameter_seed(ps))
    R=zeros(T,length(Rs),length(model.radial_spec)); dR=similar(R)
    Q=zeros(T,length(Rs),length(model.pair_spec)); dQ=similar(Q)
    rawR=similar(R); rawQ=similar(Q)
    for j in eachindex(Rs)
        jz=zs[j]
        R[j,:],dR[j,:]=radial_ed(model.radial[iz,jz],rs[j])
        Q[j,:],dQ[j,:]=radial_ed(model.pair_radial[iz,jz],rs[j])
        for n=1:size(R,2)
            rawR[j,n]=R[j,n]; scale=radial_scale(ps,:radial,n,iz,jz)
            R[j,n]*=scale; dR[j,n]*=scale
        end
        for n=1:size(Q,2)
            rawQ[j,n]=Q[j,n]; scale=radial_scale(ps,:pair,n,iz,jz)
            Q[j,n]*=scale; dQ[j,n]*=scale
        end
    end
    Y,dY=isempty(Rs) ? (zeros(T,0,(model.lmax+1)^2),zeros(SVector{3,T},0,(model.lmax+1)^2)) : angular_ed(basis,Rs)
    A=zeros(T,length(model.density))
    for (a,(n,m)) in enumerate(model.density), j in eachindex(Rs)
        A[a]+=R[j,n]*Y[j,m]
    end
    P=[prod(A[a] for a in indices) for indices in model.products]
    B=model.coupling*P
    energy=model.reference[iz]+dot(ps.WB[:,iz],B)+sum(Q[j,n]*ps.Wpair[n,iz] for j in eachindex(Rs),n=1:size(Q,2);init=zero(T))
    (;energy,Rs,rs,R,dR,Q,dQ,rawR,rawQ,Y,dY,A,P,B)
end

function site_reverse(model,ps,iz,data)
    T=eltype(data.A)
    weights=transpose(model.coupling)*ps.WB[:,iz]
    dA=zeros(T,length(model.density))
    for (k,indices) in enumerate(model.products), t in eachindex(indices)
        v=weights[k]
        for u in eachindex(indices); u==t || (v*=data.A[indices[u]]); end
        dA[indices[t]]+=v
    end
    gradients=Vector{SVector{3,T}}(undef,length(data.Rs))
    for j in eachindex(data.Rs)
        g=zero(SVector{3,T}); direction=data.Rs[j]/data.rs[j]
        for (a,(n,m)) in enumerate(model.density)
            g+=dA[a]*(data.dR[j,n]*data.Y[j,m]*direction+data.R[j,n]*data.dY[j,m])
        end
        g+=sum(ps.Wpair[n,iz]*data.dQ[j,n] for n=1:size(data.Q,2))*direction
        gradients[j]=g
    end
    (;density=dA,gradients)
end

function site_parameter_gradient(model,ps,iz,zs,data;reverse=nothing)
    T=eltype(data.A)
    gradient=map(ps) do subtree
        subtree isa NamedTuple ? map(a -> zeros(T,size(a)),subtree) : zeros(T,size(subtree))
    end
    gradient.WB[:,iz] .= data.B
    gradient.Wpair[:,iz] .= vec(sum(data.Q;dims=1))
    if hasproperty(ps,:radial)
        reverse=isnothing(reverse) ? site_reverse(model,ps,iz,data) : reverse
        for j in eachindex(zs),(a,(n,m)) in enumerate(model.density)
            gradient.radial.scale[n,iz,zs[j]]+=reverse.density[a]*data.rawR[j,n]*data.Y[j,m]
        end
    end
    if hasproperty(ps,:pair)
        for j in eachindex(zs),n=1:size(data.Q,2)
            gradient.pair.scale[n,iz,zs[j]]+=ps.Wpair[n,iz]*data.rawQ[j,n]
        end
    end
    gradient
end

function ace_vjp(args...)
    throw(ArgumentError("load ForwardDiff for local ACE higher derivative products"))
end

"""Independent serial ACE E/F/W in numeric eV/Å on fixed full environments.
Forces reverse selected products without division, so zero density factors are
valid. Every directed centre contributes once. Public outputs own their storage.
"""
function ace_evaluate(X::AbstractMatrix,model::ResolvedModel,ps,env::FullEnvironment,H::AbstractMatrix;
                      state=NamedTuple(),energy_only=false)
    size(X)==(3,length(env.species)) || throw(DimensionMismatch("positions must be 3×N"))
    size(H)==(3,3) || throw(DimensionMismatch("cell must be 3×3"))
    size(ps.WB)==(size(model.coupling,1),length(model.atomic_numbers)) || throw(DimensionMismatch("WB shape"))
    size(ps.Wpair)==(length(model.pair_spec),length(model.atomic_numbers)) || throw(DimensionMismatch("Wpair shape"))
    T=typeof(zero(eltype(X))+parameter_seed(ps))
    energy=zero(T); forces=zeros(T,size(X)); virial=zeros(T,3,3)
    basis=angular_basis(model.lmax,eltype(model.reference))
    for i in eachindex(env.species)
        data=site_diagnostics(X,model,ps,env,H,i;basis)
        energy+=data.energy
        energy_only && continue
        iz=env.species[i]
        reverse=site_reverse(model,ps,iz,data)
        for (j,k) in enumerate(data.edges)
            g=reverse.gradients[j]
            dest=env.neighbours[k]
            for a=1:3
                forces[a,i]+=g[a]; forces[a,dest]-=g[a]
                for b=1:3; virial[a,b]-=data.Rs[j][a]*g[b]; end
            end
        end
    end
    energy_only ? (;energy) : (;energy,forces,virial)
end
