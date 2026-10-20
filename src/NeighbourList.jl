"""Spatial storage and conservative cell-task schedules (no skin or cached geometry)."""
module NeighbourList
using AtomsBase, Unitful
export Workspace, prepare!, validate_schedule, PairBlock, pack

"""A 1×W callback view. Displacements are neighbour minus central; inactive geometry is (1,0,0)."""
struct PairBlock{W,V,M}
    central_species::ChemicalSpecies
    neighbour_species::NTuple{W,ChemicalSpecies}
    dx::V
    dy::V
    dz::V
    structural::M
    active::M
end

using SIMD: Vec
@inline pack(::Val{1}, f) = f(1)
@inline pack(::Val{W}, f) where W = Vec(ntuple(f, Val(W)))
@inline lane(x::Number, k) = x
@inline lane(x::Vec, k) = x[k]

"""Reusable numeric storage. Width W and Float32/Float64 precision are concrete.
`p[k]` is the original index of sorted atom k; `inverse` reverses it.
Topology and schedules are rebuilt by every `prepare!`/public evaluation.
"""
mutable struct Workspace{W,T<:Union{Float32,Float64}}
    x::Vector{T}
    y::Vector{T}
    z::Vector{T}
    species::Vector{ChemicalSpecies}
    species_map::Vector{ChemicalSpecies}
    p::Vector{Int}
    inverse::Vector{Int}
    ranges::Vector{UnitRange{Int}}
    keys::Vector{NTuple{3,Int}}
    adjacency::Vector{Vector{Int}}
    footprints::Vector{Vector{Int}}
    phases::Vector{Vector{Int}}
    lengths::NTuple{3,T}
    periodic::NTuple{3,Bool}
    cutoff::T
    force::Matrix{T}
    energies::Vector{T}
    virials::Matrix{T}
    version::UInt
    schedule_version::UInt
    timings::Dict{Symbol,Float64}
end
function Workspace(; width=4, precision=Float64)
    width in (1,2,4,8,16) || throw(ArgumentError("width must be 1, 2, 4, 8 or 16"))
    precision in (Float32,Float64) || throw(ArgumentError("precision must be Float32 or Float64"))
    _workspace(Val(width), precision)
end
function _workspace(::Val{W}, ::Type{T}) where {W,T}
    Workspace{W,T}(T[], T[], T[], ChemicalSpecies[], ChemicalSpecies[], Int[], Int[], UnitRange{Int}[],
        NTuple{3,Int}[], Vector{Int}[], Vector{Int}[], Vector{Int}[],
        (one(T),one(T),one(T)), (false,false,false), zero(T), zeros(T,3,0),
        T[], zeros(T,9,0), UInt(0), UInt(0), Dict{Symbol,Float64}())
end

# Convert units explicitly; unitless positions/cells are rejected at the public boundary.
function length_value(::Type{T}, x) where T
    x isa Unitful.Length || throw(ArgumentError("positions and cell vectors must carry length units"))
    T(ustrip(u"Å", x))
end

function prepare!(w::Workspace{W,T}, sys::AtomsBase.AbstractSystem{3}, cutoff) where {W,T}
    start = time_ns()
    c = T(cutoff)
    isfinite(c) && c > 0 || throw(ArgumentError("cutoff must be finite and positive (Å)"))
    cv = AtomsBase.cell_vectors(sys)
    H = [length_value(T, cv[j][i]) for i=1:3, j=1:3] # columns are cell vectors
    pb = Tuple(Bool.(AtomsBase.periodicity(sys)))
    for j=1:3, i=1:3
        i == j || iszero(H[i,j]) || throw(ArgumentError("only axis-aligned orthorhombic cells are supported"))
    end
    L = ntuple(i -> H[i,i], 3)
    for a=1:3
        if pb[a]
            isfinite(L[a]) && L[a] > 2c || throw(ArgumentError("periodic lengths must exceed twice cutoff"))
        end
    end
    n = length(sys)
    # Extraction is outside all callback and accumulation loops.
    coords = ntuple(_ -> Vector{T}(undef,n),3)
    species = Vector{ChemicalSpecies}(undef,n)
    empty!(w.species_map)
    for i=1:n
        r = AtomsBase.position(sys,i)
        for a=1:3
            v = length_value(T,r[a])
            isfinite(v) || throw(ArgumentError("positions must be finite"))
            coords[a][i] = pb[a] ? mod(v,L[a]) : v
        end
        sp = ChemicalSpecies(AtomsBase.species(sys,i))
        any(s -> isequal(s,sp),w.species_map) || push!(w.species_map,sp)
        species[i] = sp
    end
    w.lengths=L; w.periodic=pb; w.cutoff=c
    w.timings[:extraction] = (time_ns()-start)*1e-9
    start=time_ns()
    # Sparse occupied cells avoid allocation proportional to vacuum volume.
    origin = ntuple(a -> pb[a] || n==0 ? zero(T) : minimum(coords[a]),3)
    dims = ntuple(a -> pb[a] ? max(1,floor(Int,L[a]/c)) : 1,3)
    sizecell = ntuple(a -> pb[a] ? L[a]/dims[a] : c,3)
    atomkeys = [ntuple(a -> pb[a] ? clamp(floor(Int,coords[a][i]/sizecell[a]),0,dims[a]-1) :
        floor(Int,(coords[a][i]-origin[a])/sizecell[a]),3) for i=1:n]
    resize!(w.p,n); sortperm!(w.p,atomkeys)
    resize!(w.inverse,n)
    for v in (w.x,w.y,w.z,w.species); resize!(v,n); end
    empty!(w.keys); empty!(w.ranges)
    for k=1:n
        i=w.p[k]; w.inverse[i]=k
        w.x[k]=coords[1][i]; w.y[k]=coords[2][i]; w.z[k]=coords[3][i]; w.species[k]=species[i]
        key=atomkeys[i]
        if isempty(w.keys) || key != w.keys[end]
            push!(w.keys,key); push!(w.ranges,k:k)
        else
            w.ranges[end]=first(w.ranges[end]):k
        end
    end
    w.timings[:sorting_packing]=(time_ns()-start)*1e-9
    start=time_ns()
    lookup=Dict(key=>i for (i,key) in enumerate(w.keys))
    empty!(w.adjacency); empty!(w.footprints)
    for (i,key) in enumerate(w.keys)
        neighbours=Int[]
        for dz=-1:1, dy=-1:1, dx=-1:1
            off=(dx,dy,dz)
            other=ntuple(a -> pb[a] ? mod(key[a]+off[a],dims[a]) : key[a]+off[a],3)
            j=get(lookup,other,0)
            j>=i && push!(neighbours,j) # one unordered cell pair; same-cell i<j later
        end
        sort!(unique!(neighbours)) # periodic wraps can alias
        push!(w.adjacency,neighbours)
        push!(w.footprints,sort!(unique!([i; neighbours])))
    end
    w.version += UInt(1)
    w.timings[:topology]=(time_ns()-start)*1e-9
    start=time_ns()
    schedule!(w)
    w.timings[:schedule]=(time_ns()-start)*1e-9
    if size(w.force,2)!=n; w.force=zeros(T,3,n); end
    resize!(w.energies,length(w.keys))
    if size(w.virials,2)!=length(w.keys); w.virials=zeros(T,9,length(w.keys)); end
    w
end

function schedule!(w)
    # Cell-to-task colour incidence: no dense task×task matrix.
    incidence=[BitSet() for _ in w.keys]
    empty!(w.phases)
    for task in eachindex(w.keys)
        forbidden=BitSet()
        for cell in w.footprints[task]; union!(forbidden,incidence[cell]); end
        color=1
        while color in forbidden; color+=1; end
        if color>length(w.phases); push!(w.phases,Int[]); end
        push!(w.phases[color],task)
        for cell in w.footprints[task]; push!(incidence[cell],color); end
    end
    w.schedule_version=w.version
    w
end

"""Check phase coverage, topology version and pairwise disjoint cell write sets."""
function validate_schedule(w)
    w.schedule_version == w.version || return false
    sort(vcat(Int[],w.phases...)) == collect(eachindex(w.keys)) || return false
    for phase in w.phases
        seen=BitSet()
        for task in phase
            fp=w.footprints[task]
            isempty(intersect(seen,BitSet(fp))) || return false
            union!(seen,fp)
            all(in(fp),w.adjacency[task]) && task in fp || return false
        end
    end
    true
end
end
