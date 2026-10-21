# Initial implementation specification

## Purpose and scope

Create the reusable foundation for a Julia ACE/MACE force stack. This release
implements finite-cutoff pair potentials, spatially ordered SIMD evaluation,
energy/force/virial accumulation and conflict-free CPU cell scheduling.

The broader project aims to improve complete force-evaluation performance and
develop deep-learning implementation skills. A fast isolated pair kernel is not
the sole success criterion. Measure packing, geometry refresh, scheduling,
reduction and output permutation as well.

Existing sorting-based cell construction in
https://github.com/JuliaMolSim/NeighbourLists.jl already develops the user's idea.
Inspect its actual implementation and license. Prefer an adapter to reusable
public functionality where practical; a focused internal implementation is
acceptable if the existing representation cannot supply the needed ordering or
blocks. Document why. Do not copy substantial code without attribution and
license compliance. Do not make undocumented dependency internals a silent API.

Reference repositories:
- https://github.com/tjjarvinen/PotentialFrameWork.jl
- https://github.com/JuliaMolSim/AtomsBase.jl
- https://github.com/JuliaMolSim/AtomsSystems.jl
- https://github.com/JuliaMolSim/NeighbourLists.jl

## Suggested module boundaries

- `NeighbourList`: spatial permutation, cell ranges, candidate topology, current
  block geometry, periodic images, block traversal, write footprints and schedules.
- `PairPotentials`: callable interactions, parameter lookup, radial adapters and
  reference potential(s).
- `Execution`: traversal driver, energy/virial accumulators, endpoint force
  accumulation, output permutation and serial/parallel dispatch. It may start as
  an internal part of NeighbourList if that makes ownership clearer.

These are responsibilities, not a requirement to create empty abstraction layers.
Expose a small public API through PotentialFrameWork. Keep the root module thin.

## System boundary and storage

Accept AtomsBase-compatible systems, extracting positions, species, cell and
periodicity through the real interface. Convert once to a documented consistent
internal unit system; do not blindly strip Unitful values. Return documented
units at the public boundary. Preserve caller-owned systems and original order.

Maintain `p[k] = original index of sorted atom k` and its inverse. Sorted
coordinates, species IDs, cell ranges and force indices use the same ordering.
Use compact integer species IDs with an explicit mapping. Mapping species to
potential parameter types is owned by the potential, not guessed by the list.
Atom identity and chemical species are distinct concepts.

Store coordinates in a layout suitable for contiguous SIMD component loads.
Separate x/y/z arrays or blocked structure-of-arrays are both acceptable.
Charges and arbitrary property infrastructure are not initial deliverables.
Keep the permutation accessible so optional fields can follow it later.

## Pair blocks

Begin with one central atom and W neighbour lanes (`1 × W`). W is selectable
and concretely specialised, with W=1 as a useful reference. Do not implement
general N × M kernels before the first path works.

A block or accompanying execution metadata provides:
- One central sorted atom index and W neighbour sorted indices.
- One central species ID and W neighbour species IDs.
- Three SIMD displacement components, dx/dy/dz.
- Structural-validity mask for padding, self exclusions and absent interactions.
- An active mask accounting for physical cutoff and pair ownership.

Indices may remain executor metadata instead of callback-visible fields.
Use `r_ij = x_j + H*s_ij - x_i`. State the cell-matrix convention. Block geometry
must reflect current positions; cached displacement values must never become stale.

Retain a block when any lane is active. Evaluate masked blocks without compacting
every active pair. Padding must have valid indices and safe numerical inputs;
zeroing results after invalid division is insufficient. Never dereference sentinel
indices. Avoid accidental image/self double counting and handle diagonal blocks.

## Interaction callback and parameters

The user supplies a concrete callable object containing parameters. Conceptually:

```julia
result = interaction(block, request)
accumulate!(accumulator, metadata, block, result)
```

Names and exact dispatch are illustrative. Define and document the real API.
The callback returns per-lane energies and three per-lane Cartesian force
components. The force is on the CENTRAL endpoint. It receives enough information
to handle masking safely. The executor ensures inactive lanes contribute zero.
Support energy-only execution without calculating forces unnecessarily.

For ordinary pair potentials, enumerate each unordered physical pair once and
apply equal-and-opposite endpoint updates. If a directed list is used internally,
deduplicate or implement explicitly tested weighting/update rules. Do not mix
half-list energy accounting with full-list force accumulation.

Offer a radial adapter returning U(r) and dU/dr. With the displacement convention
above, central force is +(dU/dr)*r_ij/r. The full-block callback remains available
for specialised SIMD kernels. No dynamic dispatch per pair is required.

The potential owns parameter lookup. Support at least a mixed-species pair table
in a demonstration so the first release does not accidentally assume one species.
Compact parameter arrays and scalar central-type selection are appropriate.
General heterogeneous collections of potential objects in lanes are not.

Define physical cutoff and list radius clearly. For the first version use a
shared cutoff or a conservative list cutoff covering all interaction cutoffs.
State whether the executor or potential supplies cutoff transformations; apply
them exactly once. Supply a smooth example potential for derivative tests.

## Accumulation and virial

Evaluate and immediately consume each block result. Sum neighbour-lane forces
into the central contribution and apply negative lane forces to neighbour atoms.
Keep a central accumulator across its blocks where useful. Scatter destinations
must be valid; account for duplicate destinations before assuming vector stores
or scatters are safe. Spatial order improves locality but does not guarantee
every neighbour write is sequential.

Adopt the configurational virial convention
`W = -sum_pairs r_ij ⊗ F_i_from_j` for unordered pairs counted once.
Do not silently equate virial and stress: stress additionally involves volume and
the chosen sign convention. Derive virial in the executor from pair geometry and
forces; no nine-component per-lane result storage is necessary. Defer per-atom
virials/stresses. Energy and virial task partials can be small scalar/tensor
buffers; the one-force-array goal does not prohibit these.

Finish with one mapping back to original atom order. A permutation is sufficient;
do not sort the completed force array by a new spatial key.

## Single-array parallel cell execution

A task handles interactions owned by central atoms in one cell. Ownership rules
must prevent duplicate physical pair calculations independently of race avoidance.
Its write footprint includes every cell containing an atom it can update.

Partition tasks into phases such that footprints within a phase are disjoint.
Run a phase in parallel and synchronise before the next phase. All tasks write to
ONE shared spatially ordered force array without per-update atomics. Tasks must
also be internally race-free: keep a cell task serial over its atom/block loops
initially. Never infer safe sharing merely from assigning different central cells.

Construct and validate schedules from actual candidate-cell adjacency, including
periodic wrapping. A conservative greedy colouring of the task conflict graph
is an acceptable first implementation. Avoid constructing a dense all-cell-pair
conflict matrix: use cell-to-task incidence or other local graph construction.
Deduplicate wrapped neighbour cells. Modulo-three colouring is only an optional
fast path with a proof/check of applicability, not a universal periodic solution.
Larger interaction stencils or different cell sizes require corresponding footprints.

Tie cached schedules to the topology/version that establishes their footprints.
Do not reuse a schedule after its assumptions change. Report phase count and
parallel work per phase. Correctness is required even with only one task in a phase.
Allocate scratch per task or explicit workspace slot; do not rely on Julia task
thread IDs remaining fixed across yields/migration.

Keep task partitions cache-line aware where practical; disjoint atoms can still
false-share a cache line. This is a performance consideration, not a data race.

## Geometry and list lifetime

Initial correctness path: rebuild topology/order/schedule for each evaluation,
reusing allocations where practical. Then separate coordinate refresh from list
rebuild if implementing a skin. If a skin is supported, conservatively validate
displacement and cell changes before reuse and test periodic boundary crossings.
Do not advertise Verlet reuse without this machinery. Unsupported geometry must
raise a clear error rather than silently omit neighbours.

Minimum first-release boundary support: nonperiodic and orthorhombic periodic
systems, including mixed periodic axes. A clear minimum-image cutoff restriction
is acceptable. Triclinic and multiple-image support are extensions unless already
correctly provided by the selected backend; do not lose them accidentally if reused.

## Deferred extensions

- Optional few-force-buffer execution domains for NUMA or more concurrency.
- ACE prepared site views: active neighbours grouped by central atom and neighbour
  species, packed into uniform-species SIMD blocks. Counting/bucketing may replace
  general sorting. Preserve endpoint/image mapping for derivatives. Benchmark
  preparation overhead versus basis workload before choosing a universal policy.
- ACE site callback consumes a complete neighbourhood before producing edge
  derivatives; do not pretend the pair callback alone implements ACE.
- MACE needs whole-system layer stages and reverse feature-gradient accumulation.
  Do not assume a whole multilayer site's write footprint is one cutoff radius.
- GPU execution needs its own mapping/synchronisation design. CPU SIMD vectors
  and CPU phase scheduling are not by themselves GPU implementations.
- Long-range electrostatics, charges, training and distributed execution.

The first release must not implement these deferred features merely to fill out
the architecture. Document the extension boundaries instead.
