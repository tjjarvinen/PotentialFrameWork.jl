# PotentialFrameWork

[![Build Status](https://github.com/tjjarvinen/PotentialFrameWork.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/tjjarvinen/PotentialFrameWork.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

Finite-cutoff pair potentials for AtomsBase-compatible 3D systems, with explicit
CPU SIMD blocks and phase-scheduled cell tasks writing into one sorted force
array. Julia 1.10.10 or later is required.

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

`PairBlock` provides `central_species::Int`, `neighbour_species::NTuple{W,Int}`,
`dx/dy/dz`, `structural` and `active` masks. Species IDs index
`workspace.species_map`; extend `bind_species(p, species_map)` to prepare a concrete
parameter lookup once per evaluation. The default returns `p`. The callback must
handle all lanes; inactive geometry is replaced with `(1,0,0)` BEFORE evaluation,
and all inactive callback results are masked out by the executor. A model with a
more restricted arithmetic domain should use the masks to select suitable inputs
before division/logarithms. No invalid indices or sentinel dereferences are used.
The executor applies cutoff/ownership masks; a model's analytic taper belongs in
the potential and must be applied once. `RadialPotential(radial,cutoff)` adapts a
callable `(r, central_id, neighbour_ids, request)` returning `energy` and, for the
combined request, `derivative=dU/dr`. The adapter uses `F_i = (dU/dr) r_ij/r`.

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
AtomsCalculators integration is follow-up work: its current API was inspected,
but this initial API returns a force matrix rather than its vector of static
vectors. A future adapter can convert at the boundary. ACE site preparation,
MACE layer stages, few-buffer/NUMA domains, GPU execution, MPI, training and
long-range forces remain separate extensions.
