# PotentialFrameWork

[![Build Status](https://github.com/tjjarvinen/PotentialFrameWork.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/tjjarvinen/PotentialFrameWork.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

Finite-cutoff pair potentials for AtomsBase-compatible 3D systems, with explicit
CPU SIMD blocks and phase-scheduled cell tasks writing into one sorted force
array. Julia 1.10.10 or later and AtomsBase 0.5.3 or later are required.

```julia
using PotentialFrameWork, AtomsBase, Unitful
sys = atomic_system([:H => [0.,0.,0.]*u"Å", :He => [1.,0.,0.]*u"Å"],
    ([8.,0.,0.]*u"Å", [0.,8.,0.]*u"Å", [0.,0.,8.]*u"Å"),
    (true,true,true))
labels = species(sys, :)
potential = SmoothPair(labels, [1.0 0.7; 0.7 1.3], 2.0)
workspace = Workspace(width=4, precision=Float64)
result = evaluate!(workspace, sys, potential; parallel=true)
result.energy             # eV
result.forces             # 3 × N, original input atom order, eV/Å
result.virial             # 3 × 3, eV
energy = evaluate!(workspace, sys, potential; request=EnergyOnly()).energy
```

The pair potentials also implement the [AtomsCalculators interface](https://github.com/JuliaMolSim/AtomsCalculators.jl):

```julia
import AtomsCalculators as AC
E = AC.potential_energy(sys, potential)
F = AC.forces(sys, potential)                 # Vector{SVector{3,...}}, eV/Å
V = AC.virial(sys, potential)                 # SMatrix{3,3,...}, eV
r = AC.energy_forces_virial(sys, potential; neighbourlist=workspace, parallel=true)
f = AC.zero_forces(sys, potential)
AC.forces!(f, sys, potential; neighbourlist=workspace)  # ADDS into f
r = AC.calculate((AC.Energy(), AC.Forces(), AC.Virial()), sys, potential;
                 neighbourlist=workspace, parallel=true)
```

Omit `neighbourlist` to create a `NeighbourList.Workspace` during the call. Pass
one to reuse its storage; geometry, species, cutoff, topology and schedule are
refreshed on every call. Workspace ownership and geometry limits are the same as
for `evaluate!`. `width` and `precision` select a newly created workspace; if
specified with an existing workspace, they must match it. Default width is 4.
AtomsCalculators calls default to Float32 for Float32 system coordinates and
Float64 otherwise, including empty systems. An explicit workspace or `precision`
can override this; `zero_forces`/`zero_energy`/`zero_virial` describe the default
precision. Their units are always eV and Å.

`energy_forces`, `energy_forces_virial` and tuple `calculate` requests each use one
traversal. Energy-only calls skip force evaluation; force/virial calls use the
existing combined executor. `forces!`, `energy_forces!` and
`energy_forces_virial!` add into the caller's vector, following AtomsCalculators'
contract; clear it first when replacement is wanted. They still allocate internal
construction and boundary outputs. At this boundary forces are converted from the
existing original-order matrix into a vector of static vectors; the executor still
uses one shared sorted force array. The existing `evaluate` API keeps its matrix
output. Compatible extra keywords are accepted and ignored.

`calculate` accepts individual properties or tuples, with `parameters` and `state`
as keywords or positional arguments. Parameters remain bundled in the potential:
`nothing` or an empty parameter NamedTuple uses those values; nonempty overrides
raise an error. State is returned unchanged because the calculator is stateless.
Custom full-block potentials can subtype `AbstractPairPotential` to inherit this
interface, as shown in [the callback example](examples/custom_callback.jl).

Start Julia with `julia --threads=4 --project` to use four threads. `parallel=false`
is the default. `evaluate(sys,p; width=4, precision=Float64, ...)` creates its own
workspace. `evaluate!` reuses numeric storage; each returned result owns its output
and survives subsequent calls. A workspace must not be shared by simultaneous
public calls. Caller systems are never reordered or mutated. AtomsSystems is
optional; only AtomsBase `length`, `position(sys,i)`, `species(sys,i)`,
`cell_vectors`, and `periodicity` are required.

Coordinates/cells must carry length units and are explicitly converted to Å.
Potential parameters and callbacks use plain numeric Å and eV. `SmoothPair`
accepts actual species values (for example `species(sys,i)`), a symmetric runtime
amplitude table in eV, and cutoff in Å. It evaluates
`A_ab * (1 - r²/c²)^3` inside the strict cutoff `r<c` and zero outside. Energy and
its first two radial derivatives vanish at the cutoff. Coincident atoms are
regular for this example; singular models must define their own origin policy.
Species identity is separate from atom identity; no mixing rule is guessed.

Implementing an interaction requires no neighbour traversal or force reduction.
See the runnable [full-block and radial example](examples/custom_callback.jl).
Extend `PotentialFrameWork.cutoff(p)` with a finite positive numeric cutoff in Å,
and make `p(block, request)` callable. For `EnergyOnly()` return `(energy=U,)`;
for `EnergyForcesVirial()` return `(energy=U, fx=Fx, fy=Fy, fz=Fz)` with force on
the CENTRAL atom in eV/Å. All values are per-lane numeric values. Width 1 uses
scalar arithmetic, widths 2/4/8/16 use `SIMD.Vec`. Implement both requests to avoid
calculating forces on energy-only calls. Return values must have the geometry
precision; convert runtime parameters if needed for Float32.

`PairBlock` provides `central_species::ChemicalSpecies`,
`neighbour_species::NTuple{W,ChemicalSpecies}`,
`dx/dy/dz`, `structural` and `active` masks. Species values are preserved,
including isotope and atom-name fields. `workspace.species_map` lists the distinct
species present; extend `bind_species(p, species_map)` to prepare a concrete
parameter lookup once per evaluation. The default binding returns `p`.
`SmoothPair` uses strict `isequal` identity
for lookup (including isotope/name), rather than AtomsBase wildcard `==` matching.
Symbols passed as potential labels are converted to `ChemicalSpecies`.
The callback must handle all lanes; inactive geometry is replaced with `(1,0,0)` BEFORE evaluation,
and all inactive callback results are masked out by the executor. A model with a
more restricted arithmetic domain should use the masks to select suitable inputs
before division/logarithms. No invalid indices or sentinel dereferences are used.
The executor applies cutoff/ownership masks; a model's analytic taper belongs in
the potential and must be applied once. `RadialPotential(radial,cutoff)` adapts a
callable `(r, central_species, neighbour_species, request)` returning `energy`
and, for the combined request, `derivative=dU/dr`. The adapter uses `F_i = (dU/dr) r_ij/r`.

Cell vectors are COLUMNS of H; `r_ij = x_j + H*s_ij - x_i`. Each unordered physical
pair occurs once, with equal/opposite endpoint forces. Virial is
`W = -sum(r_ij ⊗ F_i_from_j)`; it is not stress. Under a strain ε its energy
derivative is `-sum(W .* transpose(ε))`. There is no volume normalization.

The `NeighbourList` submodule sorts atoms by cell keys into SoA coordinates and
contiguous occupied-cell ranges. `workspace.p[k]` maps sorted to original atom
index; `inverse` reverses it. Candidate adjacent cells are deduplicated after
periodic wrapping. The lower cell owns each cell pair; the lower sorted atom owns
a same-cell pair. Each serial cell task advertises all destination cells it can
update, including inactive candidates. A sparse cell-to-task greedy colouring
creates phases with disjoint footprints; each phase ends in a barrier. Forces
use one shared `3×N` array, plus small per-cell energy/virial partials. There are
no per-thread force arrays, per-update atomics, global pair-result arrays, or
thread-ID-indexed scratch buffers. Disjoint cells can still false-share cache
lines. `validate_schedule(workspace)` checks the schedule after evaluation.

Every public call rebuilds order, topology and schedule, so changed positions,
species and cells are safe. Geometry is packed from current coordinates at block
execution; there is no skin or advertised Verlet reuse. Sorted coordinate and
force buffers are reused, while construction still allocates temporary metadata.
Only axis-aligned orthorhombic cells, open axes and mixed periodic axes are
supported. Infinite diagonal lengths on open axes (isolated systems) are allowed.
Each periodic length must be STRICTLY greater than twice the cutoff; unsupported
cells raise an error. Triclinic cells and multiple-image interactions are deferred.

Run `julia --project -e 'using Pkg; Pkg.test()'`. See
[benchmark reproduction](benchmark/README.md) and the
[implementation status and measured limitations](docs/implementation/STATUS.md).
The independent all-pairs oracle, derivatives, pair ownership, actual write
footprints, widths/precisions and repeated threaded evaluations are tested.
ACE site preparation, MACE layer stages, few-buffer/NUMA domains, GPU execution,
MPI, training and long-range forces remain separate extensions.

For an external reference comparison, see the [Argon Lennard-Jones benchmark](benchmark/lennard_jones/README.md).
It matches EmpiricalPotentials' cutoff formula, checks energy/forces/virial, and
measures multiple system sizes and thread counts using AtomsSystems fixtures.
