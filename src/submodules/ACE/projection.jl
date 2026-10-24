# Truncated polynomial convolution and its division-free reverse. Only orders
# through K are stored, including intermediate orders below the retained minimum.
function projection_degree(model)
    degree=maximum(length,model.products;init=0)
    isnothing(model.policy.orders) ? degree : min(degree,last(model.policy.orders))
end
function product_projection!(P,dU,dV,U,V,products,weights,orders,fwd,rev)
    fill!(dU,0); fill!(dV,0)
    K=size(fwd,1)-1
    for (j,indices) in enumerate(products)
        fill!(fwd,0); fill!(rev,0); fwd[1,1]=1
        for (t,a) in enumerate(indices), k=0:min(t,K)
            fwd[k+1,t+1]=U[a]*fwd[k+1,t]+(k>0 ? V[a]*fwd[k,t] : zero(eltype(fwd)))
        end
        value=zero(eltype(fwd)); nu=length(indices)
        for k=0:min(nu,K)
            if isnothing(orders) || k in orders
                value+=fwd[k+1,nu+1]; rev[k+1,nu+1]=weights[j]
            end
        end
        P[j]=value
        for t=nu:-1:1
            a=indices[t]
            for k=0:min(t,K)
                g=rev[k+1,t+1]
                dU[a]+=g*fwd[k+1,t]
                rev[k+1,t]+=g*U[a]
                if k>0
                    dV[a]+=g*fwd[k,t]; rev[k,t]+=g*V[a]
                end
            end
        end
    end
    nothing
end
function projection_work(model,::Type{T}) where T
    K=projection_degree(model); nu=maximum(length,model.products;init=0)
    zeros(T,K+1,nu+1),zeros(T,K+1,nu+1)
end
function density_reverse(model,weights,U,V)
    T=promote_type(eltype(weights),eltype(U)); P=zeros(T,length(model.products))
    dU=zeros(T,length(U)); dV=similar(dU); fwd,rev=projection_work(model,T)
    product_projection!(P,dU,dV,U,V,model.products,weights,model.policy.orders,fwd,rev)
    (;P,dU,dV)
end
