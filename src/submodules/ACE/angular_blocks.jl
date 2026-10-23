# Recurrences and L2/order conventions checked against SpheriCart.jl 0.2.4,
# src/batched_kernels.jl and normalisations.jl. Adapted to explicit neighbour
# SIMD lanes with first-order Cartesian jets; no SVector{3,Vec} API assumption.
# MIT attribution: docs/implementation/licenses/SpheriCart-MIT.md.
struct AngularJet{V}
    v::V
    g::NTuple{3,V}
end
Base.:+(a::AngularJet,b::AngularJet)=AngularJet(a.v+b.v,ntuple(k -> a.g[k]+b.g[k],3))
Base.:-(a::AngularJet,b::AngularJet)=AngularJet(a.v-b.v,ntuple(k -> a.g[k]-b.g[k],3))
Base.:*(a::AngularJet,b::AngularJet)=AngularJet(a.v*b.v,ntuple(k -> a.g[k]*b.v+a.v*b.g[k],3))
Base.:*(a::Real,b::AngularJet)=AngularJet(a*b.v,ntuple(k -> a*b.g[k],3))
Base.:/(a::AngularJet,b::Real)=AngularJet(a.v/b,ntuple(k -> a.g[k]/b,3))

"""Experimental explicit SIMD neighbour blocks, with L2 real harmonics in
l*l+l+m+1 order. Temporary recurrence arrays are local per invocation, so sharing
the immutable basis across exclusive executor slots is safe. No automatic choice.
"""
struct BlockAngular{W,T}
    lmax::Int
    normalization::Matrix{T}
end
function BlockAngular(lmax; width=4,precision=Float64)
    lmax>=0 || throw(ArgumentError("lmax must be nonnegative"))
    width in (1,2,4,8,16) || throw(ArgumentError("unsupported angular width"))
    F=zeros(precision,lmax+1,lmax+1)
    for l=0:lmax
        F[l+1,1]=sqrt(precision((2l+1)/(2pi)))
        for m=1:l; F[l+1,m+1]=-F[l+1,m]/sqrt(precision((l+m)*(l+1-m))); end
    end
    BlockAngular{width,precision}(lmax,F)
end

function angular_block!(Y,dY,b::BlockAngular{W,T},Rs) where {W,T}
    n=length(Rs); L=b.lmax
    size(Y,1)>=n && size(Y,2)>=(L+1)^2 || throw(DimensionMismatch("angular outputs"))
    zero_vec=Vec{W,T}(zero(T)); one_vec=Vec{W,T}(one(T))
    jzero=AngularJet(zero_vec,(zero_vec,zero_vec,zero_vec))
    jone=AngularJet(one_vec,(zero_vec,zero_vec,zero_vec))
    Q=fill(jzero,L+1,L+1); cs=fill(jzero,L+1); sn=fill(jzero,L+1)
    for first=1:W:n
        raw=ntuple(a -> Vec(ntuple(k -> first+k-1<=n ? T(Rs[first+k-1][a]) : T(a==1),W)),3)
        radius=sqrt(raw[1]^2+raw[2]^2+raw[3]^2)
        any(radius<=0) && throw(DomainError(0,"physical zero directions are unsupported"))
        u=ntuple(a -> raw[a]/radius,3)
        jets=ntuple(a -> AngularJet(u[a],ntuple(c -> ((a==c ? one_vec : zero_vec)-u[a]*u[c])/radius,3)),3)
        x,y,z=jets; r2=x*x+y*y+z*z
        cs[1]=jone; sn[1]=jzero; Q[1,1]=jone
        for m=1:L
            cs[m+1]=cs[m]*x-sn[m]*y
            sn[m+1]=sn[m]*x+cs[m]*y
        end
        cs[1]=jone/sqrt(T(2))
        for l=0:L
            if l>0
                Q[l+1,l+1]=-(2l-1)*Q[l,l]
                Q[l+1,l]=(2l-1)*(z*Q[l,l])
                for m=0:l-2
                    Q[l+1,m+1]=((2l-1)*(z*Q[l,m+1])-(l+m-1)*(r2*Q[l-1,m+1]))/(l-m)
                end
            end
            for m=0:l
                positive=b.normalization[l+1,m+1]*(Q[l+1,m+1]*cs[m+1])
                negative=b.normalization[l+1,m+1]*(Q[l+1,m+1]*sn[m+1])
                ip=l*l+l+m+1; im=l*l+l-m+1
                for k=1:min(W,n-first+1)
                    row=first+k-1
                    Y[row,ip]=positive.v[k]
                    dY===nothing || (dY[row,ip]=SVector{3}(ntuple(a -> positive.g[a][k],3)))
                    if m>0
                        Y[row,im]=negative.v[k]
                        dY===nothing || (dY[row,im]=SVector{3}(ntuple(a -> negative.g[a][k],3)))
                    end
                end
            end
        end
    end
    Y,dY
end
angular_ed!(Y,dY,b::BlockAngular,Rs)=angular_block!(Y,dY,b,Rs)
function angular_basis_generated end
function make_angular(model,T,backend,width)
    backend===:sphericart && return angular_basis(model.lmax,T)
    backend===:generated && return angular_basis_generated(model.lmax,T)
    backend===:block && return BlockAngular(model.lmax;precision=T,width)
    throw(ArgumentError("angular_backend must be :sphericart, :generated or :block"))
end
