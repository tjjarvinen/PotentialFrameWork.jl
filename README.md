# PotentialFrameWork

[![Build Status](https://github.com/tjjarvinen/PotentialFrameWork.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/tjjarvinen/PotentialFrameWork.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

CPU pair and ACE potentials for AtomsBase-compatible 3D systems, with explicit
CPU SIMD blocks and phase-scheduled cell tasks writing into one sorted force
array. Julia 1.10.10 or later and AtomsBase 0.5.3 or later are required.

Source is grouped into three functional submodules:

```text
src/
  PotentialFrameWork.jl
  submodules/
    NeighbourList/
      NeighbourList.jl
      storage.jl
      geometry.jl
      interface.jl
      benchmarks.jl
    PairPotentials/
      PairPotentials.jl
      models.jl
      evaluation.jl
      calculators.jl
      parameters.jl
      numerical.jl
    ACE/
      ACE.jl
      resolved.jl
      environments.jl
      serial.jl
      calculators.jl
      persistence.jl
      fitting.jl
      cpu.jl
      angular_blocks.jl
```

`NeighbourList` owns spatial storage, scheduling, block geometry, neighbour queries
and neighbour-search benchmarks. `PairPotentials` owns interaction models, force
and virial evaluation, AtomsCalculators methods and parameter derivatives.
`ACE` owns resolved-model import, full species-grouped environments, independent
ACE evaluation, fitting features and CPU execution. ACEpotentials is used only
by the separate pinned reference/export scripts.
The root module reexports the public API. Low-level prepared evaluation helpers
are now in `PairPotentials`, and geometry helpers are in `NeighbourList`.

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

The geometry can also be queried through a
[NeighbourLists.jl-style interface](https://github.com/JuliaMolSim/NeighbourLists.jl):

```julia
using StaticArrays
X = [SVector(0., 0., 0.), SVector(4.9, 0., 0.)]
cell = [5. 0 0; 0 5. 0; 0 0 5.]  # lattice vectors in ROWS for this raw API
nlist = neighbour_list(X, 1.0, cell, (true, true, true))
j, R, S = neighbours(nlist, 1)
# j == [2], R[1] ≈ [-0.1, 0, 0], S[1] == [-1, 0, 0]

# AtomsBase input is supported too:
nlist = neighbour_list(sys, 2.0u"Å")
```

`j` contains original one-based indices, `R` contains Cartesian `SVector`s in
numeric Å, and `S` contains integer image `SVector`s. For raw numeric arrays,
`R[k] = X[j[k]] - X[i] + transpose(cell)*S[k]`. The row-vector convention here
matches NeighbourLists.jl; the force executor's internal H uses columns.
A 3×N position matrix is also accepted. Numeric lengths mean Å; unitful positions,
cell entries and cutoffs convert to Å. Default precision is Float64; `precision`,
`width`, `bruteforce_threshold` and `int_type=Int32` are configurable.

Each list is a materialized directed snapshot: all images within the strict
cutoff are returned, including both self-image directions, with zero-shift self
pairs excluded. Inputs and query results can be modified without changing the
snapshot; rebuild after changing geometry. Neighbour order is unspecified. This
query snapshot is separate from the `Workspace` used for force evaluation and
uses storage proportional to directed pair images. CPU construction is supported;
upstream `lazy` and GPU `backend` options are not implemented.

When NeighbourLists is loaded, an optional package extension makes
`NeighbourLists.neighbours(nlist,i)`, `NeighbourLists.nsites(nlist)` and
`NeighbourLists.npairs(nlist)` work on these snapshots too. Use qualified names
when loading both packages, for example `PotentialFrameWork.neighbour_list(...)`.
The query returns three arrays following the current NeighbourLists 0.6 API;
upstream's own types and constructors retain their original behavior. The base
interface works without installing NeighbourLists.

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
as keywords or positional arguments. `nothing` or an empty parameter NamedTuple
uses the values bundled in the potential. Explicit parameters support energy,
forces and virial; the force/virial route is a scalar correctness path. State is
returned unchanged because the pair calculator is stateless.
Custom full-block potentials can subtype `AbstractPairPotential` to inherit this
interface, as shown in [the callback example](examples/custom_callback.jl).

Energy parameter derivatives support ForwardDiff and ChainRules-compatible reverse
mode (for example Zygote). Install these packages in your application environment;
they are optional dependencies. For a `SmoothPair` potential `p`:

```julia
using ForwardDiff, Zygote, Unitful
parameters = AC.get_parameters(p)
energy(theta) = ustrip(u"eV", AC.calculate(AC.Energy(), sys, p;
    parameters=(amplitudes=theta,)).energy)
g_forward = ForwardDiff.gradient(energy, parameters.amplitudes)
g_reverse = Zygote.gradient(energy, parameters.amplitudes)[1]
```

`SmoothPair` exposes the independent upper-triangle amplitudes in column order
(`A11, A12, A22, A13, A23, A33, …`), so an off-diagonal amplitude controls both
symmetric entries. Parameter values use numeric eV/Å conventions; the example
returns derivatives of numeric energy in eV. `calculate` itself still returns
unitful energy. See the [runnable example](examples/parameter_ad.jl).

For custom potentials, implement `AC.get_parameters(p)` and
`with_parameters(p, parameters::NamedTuple)` returning a new potential. Parameters
can contain real scalars, arrays and nested NamedTuples.
Keep the physical cutoff fixed. The reconstructed energy callback must support
scalar blocks containing parameter-promoted numbers, including ForwardDiff duals.
The [benchmark LJ model](benchmark/lennard_jones/common.jl) demonstrates trainable
`epsilon` and `rmin` with its cutoff held fixed, including derivatives of its
precomputed cutoff coefficients. No implicit derivative through species, cell
construction, sorting or ownership is taken.

This parameter path uses the existing NeighbourList geometry and ownership, with
scalar energy callbacks at every workspace width. It supports `neighbourlist`,
`width`, `precision` and phased `parallel` execution, rebuilding each call; it does
not yet offer SIMD parameter derivatives. The ordinary SIMD force path is unchanged.
The numeric-eV primitive `parameter_energy(sys,p,parameters)` carries the ChainRules
`rrule`. Its gradient hook `energy_parameter_gradient` can be specialised with an
analytic implementation; loading ForwardDiff provides a generic fallback.
`SmoothPair` has an exact linear-amplitude implementation independent of
ForwardDiff, using one basis evaluation per coefficient. Reverse-mode gradients
for generic models currently use ForwardDiff locally, so their cost grows with
parameter count. Pullbacks capture gradients when the rule runs and can be reused
after workspace mutation. This legacy primitive differentiates parameters only;
its system tangent explicitly reports unsupported requests. For position and
force-output derivatives use the explicit numerical interface below. Cutoff and
cell derivatives and arbitrary higher-order reverse AD are not provided.

The explicit pair interface preserves positions and parameters through numerical
evaluation, with fixed discrete images and original atom order:

```julia
input = pair_inputs(sys, p)       # AtomsBase/unit/topology boundary
ps = AC.get_parameters(p)
result = pair_evaluate(input.positions, p, ps, input.environment, input.cell)
gx, gp = Zygote.gradient(input.positions, ps) do X, parameters
    sum(input.positions .* pair_evaluate(X, p, parameters,
        input.environment, input.cell).forces)
end
```

Both ForwardDiff and ChainRulesCore must be loaded for the supplied Zygote
pullbacks. Energy, force and virial losses differentiate into explicit positions
and parameters at fixed cell/topology. Pair force pullbacks use edge-local AD;
custom callbacks must preserve promoted scalar arithmetic. Saved pullbacks own
their inputs and survive caller/workspace mutation. Rebuild discrete metadata
after changing species, topology or cutoff. Numerical values are Å/eV throughout.

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
Orthorhombic, rotated and triclinic cells support full, mixed or no periodicity.
Cell vectors are the columns of H. Independent finite periodic vectors are
required; unused open-axis vectors may be zero or infinite. Every image within
the strict cutoff is included, even when the cutoff exceeds half a cell height.
The zero-shift self pair is excluded. Nonzero self images are counted once per
±image pair in energy and virial; their endpoint forces cancel.

`Workspace(bruteforce_threshold=32)` uses all atom pairs for systems with at most
32 atoms, avoiding spatial cell search and reducing scheduling to one task. Set
the threshold to 0 to force cell traversal for nonempty systems, or increase it to use brute force
for larger systems. The same keyword is accepted by `evaluate` and the
AtomsCalculators interface; when supplying a workspace, configure it there.
Both paths use SIMD blocks and the same periodic-image ownership. The threshold
is a configurable starting value, not a measured universal crossover. Brute force
uses one serial task even with `parallel=true`. Reciprocal-cell bounds ensure
skew-cell completeness, but highly skewed or tiny periodic cells can require
many image candidates and conservative scheduling can limit concurrency.

Compare the two strategies for your own system and cutoff:

```julia
report = benchmark_neighbourlist(system, 5.0u"Å"; parallel=true)
report.preferred          # :bruteforce or :cells
report.bruteforce.median  # seconds; raw samples are in .seconds
report.cells.median

# Force the strategy you choose for this system:
workspace = Workspace(bruteforce_threshold=length(system)) # brute force
workspace = Workspace(bruteforce_threshold=0)              # cells
```

The helper warms both paths, alternates sample order, and prints min/median/max
geometry times, preparation/traversal times and allocations. It accepts numeric
cutoffs in Å too. Options include `width`, `precision`, `parallel`, `samples=20`,
`warmup=3` and `io=nothing` to return results silently. Inputs and existing
workspaces are unchanged. Every sample includes a fresh list/schedule rebuild
and SIMD distance traversal, counting all active pair images. Potential callbacks,
force/virial accumulation and output conversion are excluded, so compare complete
calls with your real potential when choosing for production. Results apply to
the measured system and settings; run without competing tests or benchmarks.

Run `julia --project -e 'using Pkg; Pkg.test()'`. See
[benchmark reproduction](benchmark/README.md) and the
[implementation status and measured limitations](docs/implementation/STATUS.md).
The independent all-pairs oracle, derivatives, pair ownership, actual write
footprints, widths/precisions and repeated threaded evaluations are tested.
MACE layers, GPU execution, MPI and long-range forces remain separate phases.

The supplied TiAl ACE model is available through a versioned resolved export:

```julia
using SpheriCart                   # optional reference angular backend
const A = PotentialFrameWork.ACE
calc = A.ACECalculator("data/ace_ref/reference_model.json"; width=4)
result = AC.energy_forces_virial(sys, calc; parallel=true)
model, ps, st = calc.model, AC.get_parameters(calc), calc.state
input = A.ace_inputs(sys, model)
serial = A.ace_evaluate(input.positions, model, ps, input.environment, input.cell)
```

`load_recipe` checks the source hash against the resolved sibling export. To
regenerate it, run `julia --project=data/ace_ref data/ace_ref/export_resolved_model.jl`
with Julia 1.13.2. This loads the saved model without fitting. Independent runtime
evaluation has no ACEpotentials dependency. Only the resolved scalar invariant,
Agnesi/spline radial, pair and one-body components are supported. Element channels
reject isotopes and custom species names explicitly.

Load a resolved export directly with `A.ACECalculator("resolved_model_v1.toml")`.
Save a self-contained calculator checkpoint after fitting or tuning:

```julia
A.save_calculator("tuned.toml", calc)  # current installed parameters and state
# If AD training kept parameters outside the calculator:
A.save_calculator("tuned.toml", calc; parameters=trained_ps, state=st)
restored = A.ACECalculator("tuned.toml"; width=4, slots=Threads.nthreads())
```

Checkpoints include the complete resolved basis, references, current readout and
radial/pair scales, and explicit state. They need no sibling reference files or
ACEpotentials installation. Execution scratch and inference caches are rebuilt;
load with SpheriCart available for the default angular backend, or select
`angular_backend=:block`. Nested NamedTuples/tuples, numeric arrays/scalars,
strings and `nothing` are supported; Float32/Float64 types and array shapes are
preserved. Unsupported values (including AD tracers) and nonfinite numbers are
rejected. Saving replaces the destination after a complete write and does not
install an explicit parameter override in the calculator. This versioned TOML
checkpoint is a PotentialFrameWork format, separate from ACEpotentials JSON.

See the [CPU ACE interface guide](docs/ace_interface.md) for descriptor metadata,
masking/pruning, restricted construction, original-index centers, neighbor-species
order projection, custom radial callables/persistence, and optional UnicodePlots.
Runnable examples are in `examples/ace_interface.jl`.

`A.design_matrix` returns original-coordinate energy/force/virial rows and a fixed
baseline. `A.training_rows` selects present labels and residual weights. Labels
can be plain numbers (energy/virial in eV, forces in eV/Å) or Unitful quantities
in compatible units, which are converted before baseline subtraction. Force
labels accept a 3×N matrix or a vector of three-component per-atom vectors, as
used by ExtXYZ. Virial labels are configurational virials in energy units, not
volume-normalized stress. Missing properties and `nothing` omit those rows.
Positions and cells from AtomsBase systems are converted to Å by `A.ace_inputs`;
the explicit numerical kernels and design matrix use plain Å/eV arrays.
Load optional ACEfit to call `A.fit_linear`; `A.install_readout` installs the dense
original readout coordinates. `fit_pair=false` keeps imported pair weights fixed
and includes their predictions in the baseline. See the
[ACE validation and benchmark guide](benchmark/ace/README.md).

For nonlinear training, `A.trainable_parameters(model, ps)` adds explicit diagonal
radial/pair scales. Load LuxCore and call `lux_layer(model; parameters=ps)` or
`lux_layer(p)` for pair models. `LuxCore.apply(layer, input, ps, st)` returns
`(energy/forces/virial, st)` without putting mutable scratch in Lux state.
Zygote rules cover both position and parameter energy/force losses and virial
parameter losses; virial position products also work. Cell AD is explicitly
unsupported. These rules are verified for first reverse derivatives through
forces, including short fixed-cell rollouts; arbitrary reverse derivative order
and differentiation through preparation or mutable inference execution are not
claimed. Training uses the explicit numerical primitives or Lux layers.

`A.CPUExecutor` owns capacity-managed scratch slots and one sorted force array.
Rebuild the environment and call `A.prepare_executor!` before changed topology;
`A.execute_prepared!` refreshes precontraction from current parameters. Returned
`A.prepared_result` arrays own their storage. The default angular backend remains
SpheriCart. Explicit `angular_backend=:generated` / `:block` and `fused=true`
select measured experiments. Tiled fusion retains intermediates for reverse
contraction. Imported coefficients remain Float64; a Float32 executor uses
Float32 buffers with mixed radial arithmetic. No GPU support is implemented.

For an external reference comparison, see the [Argon Lennard-Jones benchmark](benchmark/lennard_jones/README.md).
It matches EmpiricalPotentials' cutoff formula, checks energy/forces/virial, and
measures multiple system sizes and thread counts using AtomsSystems fixtures.

The experimental [CPU neighborhood evaluation interface](docs/neighbour_evaluation.md)
uses local protocol declarations, explicit parameters, pair/complete-site callbacks,
and a CPU batch loop fallback. Existing pair and ACE APIs remain available;
`ACE.CPUExecutor` / `ACE.ACECalculator` accept `batch_size` (default 1).
A dependency-light [third-party callback example](examples/neighbour_evaluation.jl)
shows accessor-based pair and angular-site evaluation without a required supertype.
Pair callbacks normally consume `NeighbourInterface.distance(block)` and return
`dE/dr` with `set_radial_derivative!`; the backend supplies orientation during
force/virial reduction. Displacement access and Cartesian gradient outputs remain
available for orientation-dependent interactions.

GPU radial pair inference, device setup, ownership and limitations are documented
in [docs/gpu_pairs.md](docs/gpu_pairs.md). The default CPU APIs are unchanged.
