# Construction-time Clebsch–Gordan algebra. No coupling library or runtime
# complex arithmetic is needed by the numerical feature evaluator.
struct FeatureIrrep
    L::Int
    parity::Symbol
    channels::Int
end
function feature_irreps(spec)
    isempty(spec) && throw(ArgumentError("output_irreps must be nonempty"))
    result=map(Tuple(spec)) do s
        s isa NamedTuple && Set(keys(s))==Set((:L,:parity,:channels)) ||
            throw(ArgumentError("each irrep must be (L=..., parity=..., channels=...)"))
        s.L isa Integer && !(s.L isa Bool) && s.L>=0 || throw(ArgumentError("L must be a nonnegative integer"))
        s.channels isa Integer && !(s.channels isa Bool) && s.channels>0 || throw(ArgumentError("channels must be a positive integer"))
        s.parity in (:even,:odd) || throw(ArgumentError("parity must be :even or :odd"))
        FeatureIrrep(Int(s.L),s.parity,Int(s.channels))
    end
    length(unique((s.L,s.parity) for s in result))==length(result) || throw(ArgumentError("duplicate (L, parity) block"))
    result
end

# Condon–Shortley complex Y -> SpheriCart real L2 Y. For L=1 the
# resulting real components are sqrt(3/4π)*(y,z,x)/r.
function real_harmonic_map(l)
    U=zeros(ComplexF64,2l+1,2l+1); U[l+1,l+1]=1
    for m=1:l
        U[l+1+m,l+1+m]=(-1)^m/sqrt(2)
        U[l+1+m,l+1-m]=1/sqrt(2)
        U[l+1-m,l+1+m]=-im*(-1)^m/sqrt(2)
        U[l+1-m,l+1-m]=im/sqrt(2)
    end
    U
end
function complex_cg(a,m,b,n,L,M)
    m+n==M && abs(a-b)<=L<=a+b || return 0.0
    # BigFloat factorial arithmetic avoids integer overflow in table generation.
    f(k)=factorial(big(k))
    value=setprecision(BigFloat,256) do
        pref=sqrt(BigFloat((2L+1)*f(L+a-b)*f(L-a+b)*f(a+b-L)) / BigFloat(f(a+b+L+1)))
        pref*=sqrt(BigFloat(f(L+M)*f(L-M)*f(a-m)*f(a+m)*f(b-n)*f(b+n)))
        lo=max(0,b-L-m,a+n-L); hi=min(a+b-L,a-m,b+n)
        total=sum((BigFloat((-1)^k)/BigFloat(f(k)*f(a+b-L-k)*f(a-m-k)*f(b+n-k)*f(L-b+m+k)*f(L-a-n+k)) for k=lo:hi);init=BigFloat(0))
        pref*total
    end
    Float64(value)
end
function real_cg(a,b,L)
    Ua=real_harmonic_map(a); Ub=real_harmonic_map(b); U=real_harmonic_map(L)
    C=zeros(ComplexF64,2L+1,2a+1,2b+1)
    for m=-a:a,n=-b:b
        M=m+n; abs(M)<=L || continue
        c=complex_cg(a,m,b,n,L,M)*im^(a+b-L)
        for r=1:2L+1,s=1:2a+1,t=1:2b+1
            C[r,s,t]+=U[r,L+M+1]*c*conj(Ua[s,a+m+1])*conj(Ub[t,b+n+1])
        end
    end
    maximum(abs,imag.(C))<1e-12 || error("nonreal harmonic coupling")
    R=real.(C); R[abs.(R).<1e-14].=0
    R
end

const FeaturePolynomial=Dict{Tuple{Vararg{Int}},Float64}
function polynomial_axpy!(a,x,y)
    for (k,v) in sort!(collect(x);by=first); y[k]=get(y,k,0.0)+a*v; end
    y
end
multiplet_dot(x,y)=sum(sum(v*get(y[m],k,0.0) for (k,v) in sort!(collect(x[m]);by=first);init=0.0) for m in eachindex(x)) / length(x)
function independent_multiplet!(basis,x)
    # Reorthogonalize complete multiplets; never remove individual m components.
    for _=1:2, y in basis
        c=multiplet_dot(x,y)
        for m in eachindex(x); polynomial_axpy!(-c,y[m],x[m]); end
    end
    scale=sqrt(max(0.0,multiplet_dot(x,x)))
    scale>1e-10 || return false
    for p in x
        for k in collect(keys(p))
            v=p[k]/scale
            if abs(v)<1e-13; delete!(p,k); else; p[k]=v; end
        end
    end
    push!(basis,x); true
end
function couple_multiplet(x,y,C)
    result=[FeaturePolynomial() for _=1:size(C,1)]
    for r in axes(C,1),s in axes(C,2),t in axes(C,3)
        c=C[r,s,t]; iszero(c) && continue
        for (i,v) in sort!(collect(x[s]);by=first),(j,w) in sort!(collect(y[t]);by=first)
            key=Tuple(sort!([i...;j...])); p=result[r]
            p[key]=get(p,key,0.0)+c*v*w
        end
    end
    result
end

function generate_feature_tables(nz,rspec,correlation,lmax,totaldegree,irreps)
    # A factor is a species/radial/irrep channel. Its m components are complete.
    density=NTuple{3,Int}[]; factors=Tuple{Int,Int,Int,Int}[]; vectors=Vector{Vector{FeaturePolynomial}}()
    for z=1:nz,(r,(n,l)) in enumerate(rspec)
        0<=l<=lmax || continue
        push!(factors,(z,r,l,n+l)); indices=Int[]
        for m=-l:l
            push!(density,(z,r,l*l+l+m+1)); push!(indices,length(density))
        end
        push!(vectors,[FeaturePolynomial((i,)=>1.0) for i in indices])
    end
    blocks=[Vector{Vector{FeaturePolynomial}}() for _ in irreps]
    cache=Dict{NTuple{3,Int},Array{Float64,3}}()
    function visit(firstfactor,degree,order,parity,states)
        order>=correlation && return
        for f=firstfactor:length(factors)
            _,_,l,weight=factors[f]; degree+weight<=totaldegree || continue
            next=Dict{Int,Vector{Vector{FeaturePolynomial}}}()
            for J in sort!(collect(keys(states))),x in states[J],L=abs(J-l):J+l
                C=get!(cache,(J,l,L)) do; real_cg(J,l,L); end
                y=couple_multiplet(x,vectors[f],C)
                independent_multiplet!(get!(next,L,Vector{Vector{FeaturePolynomial}}()),y)
            end
            p=parity*(-1)^l
            for (b,s) in enumerate(irreps)
                (s.parity===:even ? 1 : -1)==p || continue
                append!(blocks[b],get(next,s.L,Vector{Vector{FeaturePolynomial}}()))
            end
            visit(f,degree+weight,order+1,p,next)
        end
    end
    visit(1,0,0,1,Dict(0=>[[FeaturePolynomial(()=>1.0)]]))
    for (b,s) in enumerate(irreps)
        isempty(blocks[b]) && throw(ArgumentError("no basis channels survive for L=$(s.L), parity=$(s.parity)"))
    end
    productkeys=Set{Tuple{Vararg{Int}}}()
    for block in blocks,x in block,p in x; union!(productkeys,keys(p)); end
    products=sort!(collect(productkeys);by=k -> (length(k),k))
    lookup=Dict(p=>i for (i,p) in enumerate(products))
    coupling=map(eachindex(irreps)) do b
        rows=Int[]; cols=Int[]; vals=Float64[]; dim=2irreps[b].L+1
        for (i,x) in enumerate(blocks[b]),m=1:dim,(k,v) in sort!(collect(x[m]);by=first)
            push!(rows,(i-1)*dim+m); push!(cols,lookup[k]); push!(vals,v)
        end
        sparse(rows,cols,vals,dim*length(blocks[b]),length(products))
    end
    density,[collect(k) for k in products],Tuple(coupling)
end
