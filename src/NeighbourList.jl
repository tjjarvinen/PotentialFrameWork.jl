"""Spatial storage and conservative cell-task schedules (no skin or cached geometry)."""
module NeighbourList
using AtomsBase, Unitful
using StaticArrays: SMatrix, SVector
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

struct CellImage{T}
    cell::Int
    shift::NTuple{3,T} # Cartesian H*s; ignored by minimum-image fast path
    self::Bool        # owns i==j only for the positive half of nonzero images
end

"""Reusable numeric storage. Width W and Float32/Float64 precision are concrete.
`p[k]` is the original index of sorted atom k; `inverse` reverses it.
Topology and schedules are rebuilt by every `prepare!`/public evaluation.
`bruteforce_threshold` (default 32 atoms) selects all-pairs traversal below or at
that size; 0 forces cells for nonempty systems. Both paths enumerate all periodic
images inside the cutoff, excluding zero-shift self pairs and owning ±self images
once. `algorithm` records `:bruteforce` or `:cells` after preparation.
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
    images::Vector{Vector{CellImage{T}}}
    minimum_image::Bool
    bruteforce_threshold::Int
    algorithm::Symbol
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
function Workspace(; width=4, precision=Float64, bruteforce_threshold=32)
    width in (1,2,4,8,16) || throw(ArgumentError("width must be 1, 2, 4, 8 or 16"))
    precision in (Float32,Float64) || throw(ArgumentError("precision must be Float32 or Float64"))
    bruteforce_threshold isa Integer && bruteforce_threshold >= 0 ||
        throw(ArgumentError("bruteforce_threshold must be a nonnegative integer"))
    _workspace(Val(width), precision, Int(bruteforce_threshold))
end
function _workspace(::Val{W}, ::Type{T}, threshold) where {W,T}
    Workspace{W,T}(T[], T[], T[], ChemicalSpecies[], ChemicalSpecies[], Int[], Int[], UnitRange{Int}[],
        NTuple{3,Int}[], Vector{Int}[], Vector{Int}[], Vector{Int}[],
        Vector{CellImage{T}}[], true, threshold, :unprepared,
        (one(T),one(T),one(T)), (false,false,false), zero(T), zeros(T,3,0),
        T[], zeros(T,9,0), UInt(0), UInt(0), Dict{Symbol,Float64}())
end

# Convert units explicitly; unitless positions/cells are rejected at the public boundary.
function length_value(::Type{T}, x) where T
    x isa Unitful.Length || throw(ArgumentError("positions and cell vectors must carry length units"))
    T(ustrip(u"Å", x))
end

# Open cell vectors do not generate images. Complete the periodic subspace with
# Cartesian directions, allowing slab/wire systems with zero or infinite vacuum
# cell vectors while requiring the actual periodic generators to be independent.
function periodic_basis(H,pb,::Type{T}) where T
    B=zeros(T,3,3)
    orthogonal=SVector{3,T}[]
    for a=1:3
        pb[a] || continue
        v=SVector{3,T}(H[:,a])
        all(isfinite,v) || throw(ArgumentError("periodic cell vectors must be finite"))
        B[:,a]=v
        for q in orthogonal; v-=sum(q.*v)*q; end
        norm=sqrt(sum(abs2,v))
        isfinite(norm) && norm>0 || throw(ArgumentError("periodic cell vectors must be linearly independent"))
        push!(orthogonal,v/norm)
    end
    for a=1:3
        pb[a] && continue
        candidates=map(1:3) do k
            v=SVector{3,T}(ntuple(j -> T(j==k),3))
            for q in orthogonal; v-=sum(q.*v)*q; end
            v
        end
        v=candidates[argmax(map(v -> sum(abs2,v),candidates))]
        v/=sqrt(sum(abs2,v))
        B[:,a]=v; push!(orthogonal,v)
    end
    SMatrix{3,3,T}(B)
end

function prepare!(w::Workspace{W,T}, sys::AtomsBase.AbstractSystem{3}, cutoff) where {W,T}
    start = time_ns()
    c = T(cutoff)
    isfinite(c) && c > 0 || throw(ArgumentError("cutoff must be finite and positive (Å)"))
    cv = AtomsBase.cell_vectors(sys)
    H = [length_value(T, cv[j][i]) for i=1:3, j=1:3] # columns are cell vectors
    pb = Tuple(Bool.(AtomsBase.periodicity(sys)))
    axisaligned = all(i==j || iszero(H[i,j]) for i=1:3,j=1:3)
    L = ntuple(i -> H[i,i], 3)
    minimum_image = !any(pb) || (axisaligned && all(!pb[a] || (isfinite(L[a]) && L[a]>2c) for a=1:3))
    # Reciprocal row norms bound every fractional component of a Cartesian ball.
    # This avoids assuming fractional rounding finds the nearest skew-cell image.
    if minimum_image
        basis = SMatrix{3,3,T}(one(T),0,0,0,one(T),0,0,0,one(T))
        reciprocal = basis
    else
        basis = periodic_basis(H,pb,T)
        reciprocal = inv(basis)
        all(isfinite,reciprocal) || throw(ArgumentError("periodic cell vectors must be linearly independent"))
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
            coords[a][i] = minimum_image && pb[a] ? mod(v,L[a]) : v
        end
        sp = ChemicalSpecies(AtomsBase.species(sys,i))
        any(s -> isequal(s,sp),w.species_map) || push!(w.species_map,sp)
        species[i] = sp
    end
    fractional = minimum_image ? coords : ntuple(_ -> Vector{T}(undef,n),3)
    if !minimum_image
        for i=1:n
            f = reciprocal*SVector(coords[1][i],coords[2][i],coords[3][i])
            f = SVector(ntuple(a -> pb[a] ? mod(f[a],one(T)) : f[a],3))
            r = basis*f
            for a=1:3; coords[a][i]=r[a]; fractional[a][i]=f[a]; end
        end
    end
    brute = n <= w.bruteforce_threshold
    w.minimum_image=minimum_image
    w.algorithm=brute ? :bruteforce : :cells
    w.lengths=L; w.periodic=pb; w.cutoff=c
    w.timings[:extraction] = (time_ns()-start)*1e-9
    start=time_ns()
    # Sparse occupied cells avoid allocation proportional to vacuum volume.
    bincoords = minimum_image ? coords : fractional
    reach = minimum_image ? ntuple(_ -> c,3) : ntuple(a -> c*sqrt(sum(abs2,reciprocal[a,:])),3)
    extent = minimum_image ? L : (one(T),one(T),one(T))
    origin = ntuple(a -> pb[a] || n==0 ? zero(T) : minimum(bincoords[a]),3)
    dims = ntuple(a -> pb[a] ? max(1,floor(Int,extent[a]/reach[a])) : 1,3)
    sizecell = ntuple(a -> pb[a] ? extent[a]/dims[a] : reach[a],3)
    atomkeys = [brute ? (0,0,0) : ntuple(a -> pb[a] ? clamp(floor(Int,bincoords[a][i]/sizecell[a]),0,dims[a]-1) :
        floor(Int,(bincoords[a][i]-origin[a])/sizecell[a]),3) for i=1:n]
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
    empty!(w.adjacency); empty!(w.footprints); empty!(w.images)
    if minimum_image
        minimum_topology!(w, dims, brute)
    else
        image_topology!(w, basis, reach, sizecell, dims, brute)
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

# The existing orthorhombic minimum-image path avoids explicit image enumeration.
function minimum_topology!(w::Workspace{W,T}, dims, brute) where {W,T}
    lookup=Dict(key=>i for (i,key) in enumerate(w.keys))
    for (i,key) in enumerate(w.keys)
        neighbours=Int[]
        if brute
            push!(neighbours,i)
        else
            for dz=-1:1, dy=-1:1, dx=-1:1
                off=(dx,dy,dz)
                other=ntuple(a -> w.periodic[a] ? mod(key[a]+off[a],dims[a]) : key[a]+off[a],3)
                j=get(lookup,other,0)
                j>=i && push!(neighbours,j)
            end
        end
        sort!(unique!(neighbours))
        push!(w.adjacency,neighbours)
        push!(w.images,[CellImage(j,(zero(T),zero(T),zero(T)),false) for j in neighbours])
        push!(w.footprints,sort!(unique!([i; neighbours])))
    end
end

# Own all shifts for i<j, and only s>0 (lexicographically) for i==j.
# Wrapped cell aliases with DIFFERENT image shifts must never be deduplicated.
function image_topology!(w::Workspace{W,T}, H, reach, sizecell, dims, brute) where {W,T}
    lookup=Dict(key=>i for (i,key) in enumerate(w.keys))
    # Wrapped fractional differences lie strictly between -1 and 1; hence
    # |s_a| < reach_a+1 and ceil(reach_a) is a sufficient integer bound.
    radius=ntuple(a -> brute ? (w.periodic[a] ? ceil(Int,reach[a]) : 0) :
        max(1,ceil(Int,reach[a]/sizecell[a])),3)
    for (i,key) in enumerate(w.keys)
        entries=CellImage{T}[]
        for z=-radius[3]:radius[3], y=-radius[2]:radius[2], x=-radius[1]:radius[1]
            off=(x,y,z)
            other=ntuple(a -> w.periodic[a] ? mod(key[a]+off[a],dims[a]) : key[a]+off[a],3)
            j=brute ? i : get(lookup,other,0)
            j>=i || continue
            shift=ntuple(a -> w.periodic[a] ? (brute ? off[a] : fld(key[a]+off[a],dims[a])) : 0,3)
            d=H*SVector(shift)
            push!(entries,CellImage(j,Tuple(d),shift>(0,0,0)))
        end
        neighbours=sort!(unique!([e.cell for e in entries]))
        push!(w.images,entries); push!(w.adjacency,neighbours)
        push!(w.footprints,sort!(unique!([i;neighbours])))
    end
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
