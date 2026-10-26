"""Independent CPU ACE numerical model and resolved-model import boundary."""
module ACE
import ..NeighbourInterface as NI
using TOML, LinearAlgebra, SparseArrays, StaticArrays
using SHA
using SIMD: Vec
using AtomsBase, Unitful
import AtomsCalculators as AC
using ..PairPotentials: force_vectors, static_virial, check_destination, add_forces!
using ..PairPotentials: parameter_seed, parameter_flatten, parameter_rebuild
using ..NeighbourList: neighbour_list, FullEnvironment, edge_displacement
import ..NeighbourList as NL

export ResolvedModel, load_resolved, radial_values, FullEnvironment, prepare_environment!, ace_inputs
export load_recipe
export ace_evaluate, site_diagnostics
export ACECalculator, design_matrix, linear_parameters, install_readout, smoothness_prior, fit_linear
export save_calculator
export training_rows
export CPUExecutor, execute_prepared!, prepare_executor!, prepared_result
export trainable_parameters, ace_vjp
export BlockAngular, angular_block!
export descriptors, basis_summary, basis_identity, BasisSelection, mask_basis, prune_basis, restrict_basis
export project_species, components, mapped_prior, coefficient_map
export restrict_parameters
export RadialMetadata, CallableRadial, replace_radial, radial_evaluate!, radial_capabilities
export allocating_radial
export register_radial!, source_radial, tabulate_radial, radial_samples, plot_basis

include("resolved.jl")
include("interface.jl")
include("radials.jl")
include("environments.jl")
include("projection.jl")
include("serial.jl")
include("calculators.jl")
include("persistence.jl")
include("interface_persistence.jl")
include("fitting.jl")
include("angular_blocks.jl")
include("cpu.jl")
end
