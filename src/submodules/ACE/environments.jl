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

function element_channel(model,sp)
    # Explicit element-only policy: reject isotopes and custom names. `==` in
    # AtomsBase uses wildcard matching, so strict identity is essential here.
    s=ChemicalSpecies(sp)
    z=Int(AtomsBase.atomic_number(s))
    isequal(s,ChemicalSpecies(z)) || throw(ArgumentError("ACE element channels do not accept isotopes or named species: $s"))
    i=findfirst(==(z),model.atomic_numbers)
    isnothing(i) && throw(ArgumentError("ACE model has no element channel for $s"))
    i
end

function prepare_environment!(env::FullEnvironment,sys,model::ResolvedModel;kwargs...)
    start=time_ns(); empty!(env.timings)
    n=length(sys); nz=length(model.atomic_numbers)
    resize!(env.species,n)
    for i=1:n; env.species[i]=element_channel(model,AtomsBase.species(sys,i)); end
    env.timings[:element_channels]=(time_ns()-start)*1e-9
    env.cutoff=maximum(ch.transform.rcut for family in (model.radial,model.pair_radial) for ch in family)
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

"""Owning AtomsBase extraction in numeric Å; H has cell vectors in columns.
Open lattice generators are zeroed before numerical evaluation.
"""
function ace_inputs(sys,model;environment=FullEnvironment(),kwargs...)
    start=time_ns()
    n=length(sys); cv=AtomsBase.cell_vectors(sys); pb=AtomsBase.periodicity(sys)
    X=reshape([ustrip(u"Å",AtomsBase.position(sys,i)[a]) for a=1:3,i=1:n],3,n)
    H=[pb[j] ? ustrip(u"Å",cv[j][a]) : zero(ustrip(u"Å",cv[j][a])) for a=1:3,j=1:3]
    elapsed=(time_ns()-start)*1e-9
    prepare_environment!(environment,sys,model;kwargs...)
    environment.timings[:boundary_extraction]=elapsed
    (;positions=X,cell=H,environment)
end

@inline function edge_displacement(X,H,env,i,k)
    # Branch before normalization/transforms; NaN times a zero mask is unsafe.
    env.active[k] || return SVector(one(eltype(X)),zero(eltype(X)),zero(eltype(X)))
    j=env.neighbours[k]
    SVector{3}(ntuple(a -> X[a,j]-X[a,i]+sum(H[a,b]*env.shifts[k][b] for b=1:3),3))
end
