"""Spatial topology, geometry, directed queries and neighbour-search diagnostics."""
module NeighbourList
import ..NeighbourInterface as NI
using AtomsBase, Unitful
import AtomsCalculators as AC
using StaticArrays: SMatrix, SVector
using SIMD: Vec, vifelse

export Workspace, PairBlock, prepare!, validate_schedule, pack,
    neighbour_list, neighbours, benchmark_neighbourlist

include("storage.jl")
include("geometry.jl")
include("interface.jl")
include("environments.jl")
include("execution.jl")
include("benchmarks.jl")
include("gpu.jl")
end
