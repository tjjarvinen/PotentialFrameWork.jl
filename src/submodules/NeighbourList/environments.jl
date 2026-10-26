"""Complete directed image graph, packed by neighbour species. All indices are
original atom indices. Padded lanes have `active=false` and safe dummy geometry.
Storage can be reused; every preparation rebuilds topology (no skin reuse).
"""
mutable struct FullEnvironment
    species::Vector{Int}
    offsets::Matrix{Int} # species boundaries per centre; padded to width
    neighbours::Vector{Int}
    shifts::Vector{SVector{3,Int}}
    active::BitVector
    footprints::Vector{Vector{Int}}
    phases::Vector{Vector{Int}}
    width::Int
    cutoff::Float64
    version::UInt
    timings::Dict{Symbol,Float64}
end
function FullEnvironment(;width=1)
    width in (1,2,4,8,16) || throw(ArgumentError("unsupported environment width"))
    FullEnvironment(Int[],zeros(Int,0,0),Int[],SVector{3,Int}[],BitVector(),Vector{Int}[],Vector{Int}[],width,0.,UInt(0),Dict{Symbol,Float64}())
end

function prepare_full_environment!(env::FullEnvironment,sys,channels,nz,cutoff;kwargs...)
    start=time_ns()
    n=length(sys)
    env.species=collect(Int,channels)
    length(env.species)==n || throw(DimensionMismatch("species mapping"))
    all(z -> 1<=z<=nz,env.species) || throw(ArgumentError("invalid species channel"))
    env.cutoff=Float64(cutoff)
    snapshot=neighbour_list(sys,env.cutoff;timings=env.timings,kwargs...)
    start=time_ns()
    size(env.offsets)==(nz+1,n) || (env.offsets=zeros(Int,nz+1,n))
    empty!(env.neighbours); empty!(env.shifts); empty!(env.active); empty!(env.footprints)
    # Counting/bucketing by species; no per-centre general sort. This baseline
    # reads the snapshot arrays directly, never the copying query API.
    for i=1:n
        js=snapshot.indices[i]; ss=snapshot.shifts[i]
        fp=Int[i]
        for z=1:nz
            env.offsets[z,i]=length(env.neighbours)+1
            count=0
            for k in eachindex(js)
                env.species[js[k]]==z || continue
                push!(env.neighbours,Int(js[k])); push!(env.shifts,SVector{3,Int}(ss[k])); push!(env.active,true)
                push!(fp,Int(js[k])); count+=1
            end
            for _=1:mod(-count,env.width)
                push!(env.neighbours,i); push!(env.shifts,SVector(0,0,0)); push!(env.active,false)
            end
        end
        env.offsets[nz+1,i]=length(env.neighbours)+1
        push!(env.footprints,sort!(unique!(fp)))
    end
    env.timings[:species_packing]=(time_ns()-start)*1e-9
    start=time_ns()
    # Full site writes touch centre AND every directed neighbour endpoint.
    # Pair lower-cell footprints are deliberately not reused.
    empty!(env.phases)
    incidence=[BitSet() for _=1:n]
    for i=1:n
        forbidden=BitSet()
        for j in env.footprints[i]; union!(forbidden,incidence[j]); end
        c=1
        while c in forbidden; c+=1; end
        c>length(env.phases) && push!(env.phases,Int[])
        push!(env.phases[c],i)
        for j in env.footprints[i]; push!(incidence[j],c); end
    end
    env.version+=UInt(1)
    env.timings[:full_environment_schedule]=(time_ns()-start)*1e-9
    env
end

@inline function edge_displacement(X,H,env,i,k)
    # Branch before normalization/transforms; NaN times a zero mask is unsafe.
    env.active[k] || return SVector(one(eltype(X)),zero(eltype(X)),zero(eltype(X)))
    j=env.neighbours[k]
    SVector{3}(ntuple(a -> X[a,j]-X[a,i]+sum(H[a,b]*env.shifts[k][b] for b=1:3),3))
end
