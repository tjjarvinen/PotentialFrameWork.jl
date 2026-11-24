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

function set_environment_species!(env,channels,n,nz,cutoff)
    length(channels)==n || throw(DimensionMismatch("species mapping"))
    resize!(env.species,n)
    for (i,z) in enumerate(channels)
        1<=z<=nz || throw(ArgumentError("invalid species channel"))
        env.species[i]=z
    end
    env.cutoff=Float64(cutoff)
end

function prepare_full_environment!(env::FullEnvironment,sys,channels,nz,cutoff;kwargs...)
    env.version+=UInt(1)
    set_environment_species!(env,channels,length(sys),nz,cutoff)
    snapshot=neighbour_list(sys,env.cutoff;timings=env.timings,kwargs...)
    pack_environment!(env,snapshot.indices,snapshot.shifts,nz,Int[])
    color_environment!(env)
    env
end

function pack_environment!(env,indices,shifts,nz,marks)
    start=time_ns()
    n=length(env.species)
    size(env.offsets)==(nz+1,n) || (env.offsets=Matrix{Int}(undef,nz+1,n))
    empty!(env.neighbours); empty!(env.shifts); empty!(env.active)
    refill_vectors!(env.footprints,n)
    resize!(marks,n); fill!(marks,0)
    # Preserve the traversal order within each species bucket.
    for i=1:n
        js=indices[i]; ss=shifts[i]
        fp=env.footprints[i]
        push!(fp,i); marks[i]=i
        for z=1:nz
            env.offsets[z,i]=length(env.neighbours)+1
            count=0
            for k in eachindex(js)
                j=js[k]
                env.species[j]==z || continue
                push!(env.neighbours,j); push!(env.shifts,ss[k]); push!(env.active,true)
                if marks[j]!=i
                    marks[j]=i
                    push!(fp,j)
                end
                count+=1
            end
            for _=1:mod(-count,env.width)
                push!(env.neighbours,i); push!(env.shifts,SVector(0,0,0)); push!(env.active,false)
            end
        end
        env.offsets[nz+1,i]=length(env.neighbours)+1
        sort!(fp)
    end
    env.timings[:species_packing]=(time_ns()-start)*1e-9
    env
end

function color_environment!(env)
    start=time_ns()
    # Full site writes touch centre AND every directed neighbour endpoint.
    # Pair lower-cell footprints are deliberately not reused.
    empty!(env.phases)
    incidence=[BitSet() for _ in env.species]
    for i in eachindex(env.species)
        forbidden=BitSet()
        for j in env.footprints[i]; union!(forbidden,incidence[j]); end
        c=1
        while c in forbidden; c+=1; end
        c>length(env.phases) && push!(env.phases,Int[])
        push!(env.phases[c],i)
        for j in env.footprints[i]; push!(incidence[j],c); end
    end
    env.timings[:full_environment_schedule]=(time_ns()-start)*1e-9
    env
end

@inline function edge_displacement(X,H,env,i,k)
    # Branch before normalization/transforms; NaN times a zero mask is unsafe.
    env.active[k] || return SVector(one(eltype(X)),zero(eltype(X)),zero(eltype(X)))
    j=env.neighbours[k]
    SVector{3}(ntuple(a -> X[a,j]-X[a,i]+sum(H[a,b]*env.shifts[k][b] for b=1:3),3))
end
