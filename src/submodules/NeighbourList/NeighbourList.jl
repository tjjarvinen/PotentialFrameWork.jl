"""Spatial topology, geometry, directed queries and neighbour-search diagnostics."""
module NeighbourList
using AtomsBase, Unitful
using StaticArrays: SMatrix, SVector
using SIMD: Vec, vifelse

export Workspace, PairBlock, prepare!, validate_schedule, pack,
    neighbour_list, neighbours, benchmark_neighbourlist

include("storage.jl")
include("geometry.jl")
include("interface.jl")
include("benchmarks.jl")
end
