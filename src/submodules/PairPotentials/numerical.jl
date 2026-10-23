"""Fixed unordered pair-image metadata in original atom order, without saved
displacements. Build again after topology, species or cutoff changes. `cell` is
passed separately to `pair_evaluate` with lattice vectors in columns.
"""
struct PairEnvironment{R}
    labels::Vector{ChemicalSpecies}
    centres::Vector{Int}
    neighbours::Vector{Int}
    shifts::Vector{SVector{3,Int}}
    radius::R
end

"""Extract owning numeric Å positions/cell and fixed pair metadata from AtomsBase.
Topology construction is discrete; differentiate `pair_evaluate`, with its explicit
positions and parameters, after this boundary. No workspace is retained.
"""
function pair_inputs(sys::AbstractSystem{3}, p; kwargs...)
    X=reshape([ustrip(u"Å",position(sys,i)[a]) for a=1:3,i=1:length(sys)],3,length(sys))
    cv=AtomsBase.cell_vectors(sys)
    H=[ustrip(u"Å",cv[j][a]) for a=1:3,j=1:3]
    # Open cell vectors may contain Inf; their zero shifts must never multiply it.
    pb=AtomsBase.periodicity(sys)
    H=[pb[j] ? H[a,j] : zero(H[a,j]) for a=1:3,j=1:3]
    snapshot=neighbour_list(sys,cutoff(p);kwargs...)
    centres=Int[]; dests=Int[]; shifts=SVector{3,Int}[]
    for i=1:length(sys), k in eachindex(snapshot.indices[i])
        j=Int(snapshot.indices[i][k]); s=SVector{3,Int}(snapshot.shifts[i][k])
        (i<j || (i==j && Tuple(s)>(0,0,0))) || continue
        push!(centres,i); push!(dests,j); push!(shifts,s)
    end
    env=PairEnvironment(ChemicalSpecies[AtomsBase.species(sys,i) for i=1:length(sys)],
        centres,dests,shifts,cutoff(p))
    (;positions=X,cell=H,environment=env)
end

# Scalar explicit arithmetic must preserve both parameter and position types.
# Ordinary SIMD inference retains its existing geometry-precision policy.
numerical_pair(p,b,request)=p(b,request)
function numerical_pair(p::BoundSmooth,b,request)
    a=p.amplitudes[p.indices[b.central_species],p.indices[b.neighbour_species[1]]]
    q=1-(b.dx*b.dx+b.dy*b.dy+b.dz*b.dz)/p.radius^2
    energy=a*q^3
    request isa EnergyOnly && return (;energy)
    g=-6a*q^2/p.radius^2
    (;energy,fx=g*b.dx,fy=g*b.dy,fz=g*b.dz)
end

function edge_result(bound,env,k,d,request)
    b=PairBlock(env.labels[env.centres[k]],(env.labels[env.neighbours[k]],),
        d[1],d[2],d[3],true,true)
    numerical_pair(bound,b,request)
end

"""Numeric eV/Å pair E/F/W with explicit `3×N` positions, parameters and `3×3`
column cell. Environment indices/images are fixed, geometry is recomputed.
Returns owning outputs. Scalar serial correctness path; one force array, no
per-worker forces. Requires a valid environment for the current positions.
`state` contains no mutable numerical scratch. Cell differentiation is unsupported.
"""
function pair_evaluate(X::AbstractMatrix,p::AbstractPairPotential,ps::NamedTuple,
                       env::PairEnvironment,H::AbstractMatrix;request=EnergyForcesVirial(),state=NamedTuple())
    size(X)==(3,length(env.labels)) || throw(DimensionMismatch("positions must be 3×N"))
    size(H)==(3,3) || throw(DimensionMismatch("cell must be 3×3"))
    request isa Union{EnergyOnly,EnergyForcesVirial} || throw(ArgumentError("unsupported request"))
    updated=with_parameters(p,ps)
    cutoff(updated)==env.radius==cutoff(p) || throw(ArgumentError("environment cutoff mismatch"))
    bound=bind_species(updated,unique(env.labels))
    seed=zero(eltype(X))+zero(eltype(H))+parameter_seed(ps)
    energy=seed
    forces=fill(seed, size(X)); virial=fill(seed,3,3)
    for k in eachindex(env.centres)
        i=env.centres[k]; j=env.neighbours[k]
        d=SVector{3}(ntuple(a -> X[a,j]-X[a,i]+sum(H[a,b]*env.shifts[k][b] for b=1:3)+seed,3))
        sum(abs2,d)<env.radius^2 || continue
        out=edge_result(bound,env,k,d,request)
        energy+=out.energy
        if request isa EnergyForcesVirial
            g=(out.fx,out.fy,out.fz)
            for a=1:3
                forces[a,i]+=g[a]; forces[a,j]-=g[a]
                for b=1:3; virial[a,b]-=d[a]*g[b]; end
            end
        end
    end
    request isa EnergyOnly ? (;energy) : (;energy,forces,virial)
end

"""Local derivative-product hook, supplied by the ForwardDiff extension.
This correctness fallback computes edge-local products, never a full Hessian.
"""
function pair_vjp(args...)
    throw(ArgumentError("load ForwardDiff to compute pair derivative products"))
end

"""Create a LuxCore layer when LuxCore is loaded. Inputs hold explicit positions,
cell and fixed environment; state never contains reusable mutable scratch.
"""
function lux_layer(args...;kwargs...)
    throw(ArgumentError("load LuxCore to create a potential layer"))
end
