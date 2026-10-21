"""Materialized directed neighbour snapshot in original atom order.
Use `neighbour_list` and `neighbours`; queries return independent arrays.
"""
struct NeighbourSnapshot{T,TI}
    indices::Vector{Vector{TI}}
    displacements::Vector{Vector{SVector{3,T}}}
    shifts::Vector{Vector{SVector{3,TI}}}
end

# Keep raw-array handling at the AtomsBase boundary rather than in hot loops.
struct CoordinateSystem{T} <: AtomsBase.AbstractSystem{3}
    positions::Vector{SVector{3,T}}
    H::SMatrix{3,3,T,9}
    pbc::NTuple{3,Bool}
end
Base.length(s::CoordinateSystem)=length(s.positions)
AtomsBase.position(s::CoordinateSystem,i::Integer)=s.positions[i]*u"Å"
AtomsBase.cell_vectors(s::CoordinateSystem)=ntuple(a -> s.H[:,a]*u"Å",3)
AtomsBase.periodicity(s::CoordinateSystem)=s.pbc
AtomsBase.species(::CoordinateSystem,::Integer)=ChemicalSpecies(:H)

function numeric_length(::Type{T},v) where T
    v isa Unitful.Length && return T(ustrip(u"Å",v))
    v isa Real || throw(ArgumentError("lengths must be real numbers or unitful lengths"))
    T(v)
end

"""
    neighbour_list(X, cutoff, cell, pbc; width=4, precision=Float64,
                   bruteforce_threshold=32, int_type=Int32)
    neighbour_list(sys::AtomsBase.AbstractSystem{3}, cutoff; kwargs...)

Build a materialized directed neighbour snapshot. `X` is a vector of 3-vectors
(or a 3×N matrix); `cell` is a 3×3 matrix with lattice vectors in ROWS, matching
NeighbourLists.jl's raw-array API. Numeric lengths are in Å; unitful lengths are
converted to Å. For AtomsBase input, cell vectors follow its usual convention.
Inputs are copied; rebuild the snapshot after changing positions/cell/cutoff.

`neighbours(list,i)` returns `(js,Rs,Ss)`: original one-based atom indices,
Cartesian `SVector` displacements in numeric Å, and integer `SVector` shifts.
For raw input, `Rs[k] = X[js[k]] - X[i] + transpose(cell)*Ss[k]` (in Å).
All images strictly inside the cutoff are included, including both ±self images;
only the zero-shift self pair is excluded. Neighbour ordering is unspecified.
The index/shift type defaults to Int32; overflow raises an error. These CPU
snapshots do not implement upstream `lazy` or GPU `backend` options.
"""
function neighbour_list(X::Union{AbstractVector,AbstractMatrix},cutoff,
                        cell::AbstractMatrix,pbc;
                        width=4,precision=Float64,bruteforce_threshold=32,int_type=Int32)
    size(cell)==(3,3) || throw(DimensionMismatch("cell must be 3×3 with vectors in rows"))
    length(pbc)==3 || throw(DimensionMismatch("pbc must contain three flags"))
    X isa AbstractMatrix && size(X,1)!=3 && throw(DimensionMismatch("positions must be a 3×N matrix"))
    n=X isa AbstractMatrix ? size(X,2) : length(X)
    # Validate options before numeric conversion, using the same workspace policy.
    Workspace(;width,precision,bruteforce_threshold)
    positions=Vector{SVector{3,precision}}(undef,n)
    for i=1:n
        r=X isa AbstractMatrix ? view(X,:,i) : X[i]
        length(r)==3 || throw(DimensionMismatch("each position must have three components"))
        positions[i]=SVector{3,precision}(ntuple(a -> numeric_length(precision,r[a]),3))
    end
    H=SMatrix{3,3,precision}(ntuple(k -> numeric_length(precision,cell[cld(k,3),mod1(k,3)]),9))
    sys=CoordinateSystem(positions,H,Tuple(Bool.(pbc)))
    neighbour_list(sys,cutoff;width,precision,bruteforce_threshold,int_type)
end

function neighbour_list(sys::AtomsBase.AbstractSystem{3},cutoff;
                        width=4,precision=Float64,bruteforce_threshold=32,int_type=Int32)
    int_type in (Int8,Int16,Int32,Int64,Int128) ||
        throw(ArgumentError("int_type must be a signed built-in integer type"))
    n=length(sys)
    n<=typemax(int_type) || throw(ArgumentError("int_type cannot represent all atom indices"))
    w=Workspace(;width,precision,bruteforce_threshold)
    c=numeric_length(precision,cutoff)
    prepare!(w,sys,c)
    materialize(w,sys,int_type)
end

function materialize(w::Workspace{W,T},sys,::Type{TI}) where {W,T,TI}
    n=length(sys)
    indices=[TI[] for _=1:n]
    displacements=[SVector{3,T}[] for _=1:n]
    shifts=[SVector{3,TI}[] for _=1:n]
    # Undo the canonical wrapping performed by prepare!, retaining exact image
    # labels relative to caller positions even several cells outside the box.
    pb=w.periodic
    cv=AtomsBase.cell_vectors(sys)
    H=SMatrix{3,3,T}([numeric_length(T,cv[j][a]) for a=1:3,j=1:3])
    reciprocal=w.minimum_image ? nothing : inv(periodic_basis(H,pb,T))
    winding=Vector{SVector{3,Int}}(undef,n)
    for i=1:n
        position=AtomsBase.position(sys,i)
        r=SVector{3,T}(ntuple(a -> numeric_length(T,position[a]),3))
        frac=w.minimum_image ? r : reciprocal*r
        winding[i]=SVector{3,Int}(ntuple(a -> !pb[a] ? 0 :
            floor(Int,w.minimum_image ? frac[a]/w.lengths[a] : frac[a]),3))
    end
    for task in eachindex(w.ranges), i in w.ranges[task], image in w.images[task]
        r=w.ranges[image.cell]
        for j=max(first(r),i+(image.self ? 0 : 1)):W:last(r)
            js,b=geometry(w,i,j,last(r),image)
            for k=1:W
                lane(b.active,k) || continue
                dest=js[k]; a=w.p[i]; bb=w.p[dest]
                R=SVector(lane(b.dx,k),lane(b.dy,k),lane(b.dz,k))
                lattice=if w.minimum_image
                    d=(w.x[dest]-w.x[i],w.y[dest]-w.y[i],w.z[dest]-w.z[i])
                    SVector{3,Int}(ntuple(axis -> pb[axis] ? -round(Int,d[axis]/w.lengths[axis]) : 0,3))
                else
                    SVector(image.lattice)
                end
                S=SVector{3,Int}(ntuple(axis -> Base.Checked.checked_sub(
                    Base.Checked.checked_add(lattice[axis],winding[a][axis]),winding[bb][axis]),3))
                shift=SVector{3,TI}(S) # checked conversion; no silent narrowing
                opposite=SVector{3,TI}(ntuple(axis ->
                    Base.Checked.checked_sub(zero(TI),shift[axis]),3))
                push!(indices[a],TI(bb)); push!(displacements[a],R); push!(shifts[a],shift)
                push!(indices[bb],TI(a)); push!(displacements[bb],-R); push!(shifts[bb],opposite)
            end
        end
    end
    NeighbourSnapshot(indices,displacements,shifts)
end

"""Return `(js,Rs,Ss)` for original atom `i`, as independent arrays.
Displacements are numeric Å; shifts are dimensionless integers. Bounds are checked.
"""
function neighbours(list::NeighbourSnapshot,i::Integer)
    checkbounds(list.indices,i)
    copy(list.indices[i]),copy(list.displacements[i]),copy(list.shifts[i])
end
