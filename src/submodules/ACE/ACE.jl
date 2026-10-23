"""Independent CPU ACE numerical model and resolved-model import boundary."""
module ACE
using TOML, LinearAlgebra, SparseArrays, StaticArrays
using SHA
using SIMD: Vec
using AtomsBase, Unitful
import AtomsCalculators as AC
using ..PairPotentials: force_vectors, static_virial, check_destination, add_forces!
using ..PairPotentials: parameter_seed, parameter_flatten, parameter_rebuild
using ..NeighbourList: neighbour_list

export ResolvedModel, load_resolved, radial_values, FullEnvironment, prepare_environment!, ace_inputs
export load_recipe
export ace_evaluate, site_diagnostics
export ACECalculator, design_matrix, linear_parameters, install_readout, smoothness_prior, fit_linear
export save_calculator
export training_rows
export CPUExecutor, execute_prepared!, prepare_executor!, prepared_result
export trainable_parameters, ace_vjp
export BlockAngular, angular_block!

include("resolved.jl")
include("environments.jl")
include("serial.jl")
include("calculators.jl")
include("persistence.jl")
include("fitting.jl")
include("angular_blocks.jl")
include("cpu.jl")
end
