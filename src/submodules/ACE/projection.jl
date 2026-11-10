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
# Prepared scalar inference only. Functional AD and feature evaluation retain
# product_projection! above. The first row of existing bounded scratch stores
# prefixes; reverse needs one scalar, including at zero/repeated factors.
function unrestricted_products!(P,dU,dV,U,V,products,weights,prefix)
    fill!(dU,0); fill!(dV,0)
    for (j,indices) in enumerate(products)
        prefix[1,1]=one(eltype(prefix))
        for (t,a) in enumerate(indices)
            prefix[1,t+1]=prefix[1,t]*(U[a]+V[a])
        end
        P[j]=prefix[1,length(indices)+1]
        adj=weights[j]
        for t=length(indices):-1:1
            a=indices[t]
            g=adj*prefix[1,t]
            dU[a]+=g; dV[a]+=g
            adj*=U[a]+V[a]
        end
    end
    nothing
end
function check_product_algorithm(algorithm)
    algorithm in (:auto,:reference,:optimized) ||
        throw(ArgumentError("product_algorithm must be :auto, :reference or :optimized"))
    algorithm
end
function prepared_products!(P,dU,dV,U,V,products,weights,orders,fwd,rev,algorithm)
    if algorithm!==:reference && isnothing(orders)
        unrestricted_products!(P,dU,dV,U,V,products,weights,fwd)
    else
        product_projection!(P,dU,dV,U,V,products,weights,orders,fwd,rev)
    end
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
