# Explicit numerical evidence for canonical factor reordering with nonzero U/V.
include("common.jl")
function main()
    products=[[2],[2,1],[1,2],[2,2],[2,1,3],[4,3,2,1],
              [4,3,2,1,4],[2,1,2,3,4,4],[1,1],[1,1,1]]
    plan=A.ProductDAG(products,5); records=Dict[]
    for T in (Float32,Float64), algorithm in (:reference,:optimized)
        tol=T===Float32 ? 2e-5 : 3e-12
        maxima=zeros(3); scaled=zeros(3); calls=0
        for u in ([.7,-.3,1.2,-.6,3.],[0.,-.3,1.2,.6,0.],[0.,0.,1.2,-.6,2.]),
            v in (zeros(5),[.2,.3,-.4,.1,0.]), pass in 1:3
            U=T.(u).*T(pass); V=T.(v)
            weights=T[isodd(j+pass) ? 0 : sin(j)/pass for j in eachindex(products)]
            P=zeros(T,length(products)); dU=zeros(T,5); dV=similar(dU)
            fwd=zeros(T,7,7); rev=similar(fwd)
            A.prepared_products!(P,dU,dV,U,V,products,weights,nothing,fwd,rev,algorithm)
            reference=(copy(P),copy(dU),copy(dV))
            values=zeros(T,5+length(plan.left)); adj=similar(values)
            A.dag_products!(P,dU,dV,U,V,weights,plan,values,adj)
            for (i,(y,r)) in enumerate(zip((P,dU,dV),reference))
                @assert isapprox(y,r;atol=tol,rtol=tol)
                delta=abs.(Float64.(y)-Float64.(r))
                maxima[i]=max(maxima[i],maximum(delta))
                scaled[i]=max(scaled[i],maximum(delta./(tol .+ tol.*abs.(Float64.(r)))))
            end
            calls+=1
        end
        push!(records,Dict("precision"=>string(T),"comparison"=>string(algorithm),"calls"=>calls,
            "P_dU_dV_max_abs"=>maxima,"P_dU_dV_max_scaled_element_error"=>scaled,"atol"=>tol,"rtol"=>tol))
    end
    open(io->TOML.print(io,Dict("julia"=>string(VERSION),"script_sha256"=>bytes2hex(sha256(read(@__FILE__))),"cases"=>records);sorted=true),ARGS[1],"w")
end
main()
