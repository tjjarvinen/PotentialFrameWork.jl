"""Spatial topology, geometry, directed queries and neighbour-search diagnostics."""
module NeighbourList
import ..NeighbourInterface as NI
using AtomsBase, Unitful
import AtomsCalculators as AC
using StaticArrays: SMatrix, SVector
using SIMD: Vec, vifelse

export ChunkBuffers, prepare_neighbors, form_numa_workers, foreach_worker!, get_buffers, get_scratch,
    get_results, export_field, refresh_positions!, ChunkTimings, site_ordinal, central_slot, neighbor_slot, central_feature, neighbor_feature

export Workspace, PairBlock, prepare!, validate_schedule, pack,
    neighbour_list, neighbours, benchmark_neighbourlist

include("storage.jl")
include("geometry.jl")
include("interface.jl")
include("environments.jl")
include("execution.jl")
include("chunk_threads.jl")
include("chunks.jl")
include("pullbacks.jl")
export site_evaluate, site_pullback
include("benchmarks.jl")
include("gpu.jl")
end
