"""Experimental local neighborhood protocol. No backend or model dependency.

Callbacks overwrite requested owned output slots and return `nothing`. Pair
callbacks normally return dE/dr from distance-only mathematics. Optional Cartesian
pair gradients and site gradients are ∂E/∂(Xj-Xi+H*s), never global atom forces. Views are borrowed through callback
and immediate backend consumption and must not be retained. Inactive geometry
must be safe before mathematical operations; inactive outputs are not consumed.
Properties are canonical tuples of AtomsCalculators property instances. Numerical
units are declared with AtomsCalculators.length_unit/energy_unit on potentials.
See docs/neighbour_evaluation.md for preparation and accessor contracts.
"""
module NeighbourInterface
import AtomsCalculators
export evaluate_pair!, evaluate_site!, evaluate_sites!, evaluation_kind, search_cutoff,
    species_channels, prepare_potential, update_prepared!, make_scratch, central_index, central_species,
    neighbour_count, neighbour_species, distance, displacement, active, structural, species_range,
    site_count, site_at, output_at, set_energy!, set_radial_derivative!, set_gradient!, output_energy, output_radial_derivative, output_gradient

"""Overwrite requested pair energies and radial derivatives dE/dr (the usual path).
Use distance(block) for scalar/SIMD distances. Optional displacement(block) and
Cartesian set_gradient! support orientation-dependent pairs. Backend owns radial
orientation conversion and global reduction; callbacks return nothing."""
function evaluate_pair! end
"""Overwrite a complete center's requested energy and directed-edge gradients."""
function evaluate_site! end
"""Overwrite independent outputs for complete sites; backend may supply a loop fallback."""
function evaluate_sites! end
"""Return :pair or :site. Resolved once at preparation, not per interaction."""
function evaluation_kind end
"""Physical search radius in the potential's declared numeric length unit."""
function search_cutoff end
"""Ordered required species channels, or nothing for provider-native species."""
function species_channels end
"""prepare_potential(p, species; precision, width, kwargs...) -> concrete evaluator.
Validate layout, species mapping and units here; caller parameters stay authoritative."""
function prepare_potential end
"""make_scratch(prepared, capacity, precision) -> exclusive reusable local scratch."""
function make_scratch end
"""Refresh optional derived coefficients from current explicit parameters before traversal."""
function update_prepared! end
function central_index end
function central_species end
function neighbour_count end
function neighbour_species end
"""displacement(view, edge) -> numeric three-vector; pair form displacement(block) -> components."""
function displacement end
"""Pair distance in declared length units (scalar or independent SIMD lanes).
This is the primary pair geometry accessor; displacement is optional for radial
callbacks. Inactive distances are safe before evaluation. At a physical origin,
radial models must have dE/dr=0 or explicitly reject that geometry."""
function distance end
function active end
function structural end
"""Local edge range for a required channel, including explicitly inactive padding."""
function species_range end
function site_count end
function site_at end
function output_at end
function set_energy! end
function set_gradient! end
"""Set pair dE/dr, not force magnitude -dE/dr. Backend computes (dE/dr)*d/r.
Use this for radial pairs; set_gradient! remains the optional Cartesian path."""
function set_radial_derivative! end
function output_energy end
function output_gradient end
function output_radial_derivative end
end
