# GPU phase 1 specification

## 1. Scope and boundaries

Deliver a portable GPU sorted-cell neighbor backend and generic radial pair
executor, primarily validated on AMDGPU with an RX 6800 XT. GPU ACE/MACE and
GPU differentiation are later phases. Existing CPU behavior, periodic images,
functional AD, model persistence and public APIs remain supported.

Use KernelAbstractions (KA) for custom kernels and AcceleratedKernels (AK) for
sorting/scans/reductions where appropriate. Verify the installed versions' APIs;
do not copy signatures from unversioned documentation without checking them.
Do not introduce mandatory AMDGPU/CUDA loading or require a GPU to load the core
package. Use a coherent optional extension/dependency design based on current
Julia extension and backend APIs. No bespoke backend-specific force kernel is
required. Keep the core Julia compatibility where practical.

Current inspected boundaries:
- `src/submodules/NeighbourInterface/NeighbourInterface.jl`: pure declarations.
- `NeighbourList/storage.jl`: CPU Workspace, sorted SoA, permutations, cell/image
  topology and coloring; CPU fields include host-only nested vectors/dictionaries.
- `NeighbourList/execution.jl`: callback views, property normalization, reduction.
- `PairPotentials/models.jl`: radial callbacks, prepared adapters and species binding.
- `PairPotentials/evaluation.jl`: public preparation/execution/output boundary.
- `PairPotentials/numerical.jl` and extensions: existing functional CPU AD.

Reinspect these rather than assuming the local checkout is unchanged. Add a
separate GPU workspace/execution representation; do not force CuArray/ROCArray
into fields concretely declared as CPU Vector/Matrix. Keep the root module thin.

## 2. Configuration and results

Integrate backend selection and execution strategy with the completed interface.
Choose exact keyword/type names after inspection; document executable examples.
Provide explicit Float32/Float64 computation precision independent of system
precision. Preserve defaults on existing CPU APIs. For the new GPU entry point,
document its default; preserving input Float32/Float64 is preferred when unspecified.

Strategies: central ownership (default), candidate-thread tiling (experimental).
Skin defaults to zero/disabled. Strategy, capacity and launch dimensions are
execution configuration, not physical potential parameters. CPU SIMD width is
not a GPU workgroup size. Make launch/tile dimensions independently tunable.

Keep an owning device-numerical result path for repeated execution and benchmarking.
The existing public convenience result keeps its units, requested property keys,
owning outputs and original atom order. Document where host copies/synchronization
occur. No result may unexpectedly alias reusable workspace after a later call.
Do not copy forces to the CPU in a device-resident execution solely to satisfy
the convenience API. Respect all nonempty Energy/Forces/Virial subsets and the
existing empty-request semantics.

## 3. Precision, species and authoritative parameters

Selected computation type governs positions, cell, cutoff/skin, reference
positions, floating-point parameter tables, intermediates, energy, forces and
virial. Initially accumulation uses that same type. Integer indexing is separate;
use checked bounds and do not silently overflow atom/cell/entry offsets. Avoid
Float64 literals promoting Float32 kernels. Preparation never mutates the system.

Retain full ChemicalSpecies identity. The user designed its 64-bit representation
for integer conversion: verify the actual dependency's supported conversion,
round-trip, equality and device-compatible operations. Do not equate species with
atomic number or discard named/isotopic identity. Device records may carry the
verified fixed-width identity representation. Public accessors must preserve the
existing species semantics rather than silently returning unrelated channel IDs.

Existing dictionary-based potential setup stays valid on the host. Preparation
may additionally produce dense integer channels and device parameter tables.
No ordinary Julia Dict or string lookup in kernels. Test named/isotopic identity
where supported by the potential; reject unsupported model species explicitly.

Every call's explicit parameters remain authoritative. Preparing tables once must
not make subsequent parameter updates stale. Provide documented refresh behavior
and invalidation on model species/cutoff changes. Geometry precision conversion
here is inference behavior, not an AD implementation. Preserve independent CPU
parameter precision/AD behavior unchanged.

## 4. Sorted storage and cell construction

Retain the original record-sort approach: records carry position, species identity,
original atom index and cell key. Sort by cell key with AK; then unpack into aligned
sorted coordinate/species/index arrays used directly by kernels. Use existing SoA
coordinate organization initially if suitable. Preserve an inverse permutation
where needed. Sorting key-only and gathering is not a required optimization: the
user measured sorting as negligible. Measure sorting, do not presume it dominates.

Use cell keys derived from actual dimensions with checked representability; no
fixed base-1000 packing. Handle empty cells and empty systems explicitly. Current
CPU storage is sparse in occupied cells; avoid a dense allocation proportional to
arbitrarily large vacuum. Prefer sorted unique occupied keys plus ranges and a
device-compatible lookup, or document a bounded dense mode with sparse fallback.
Cell boundaries can be detected from sorted-key transitions; AK compaction/scan
may help. An atomic histogram is an alternative only if measurements justify it.
Keep cell construction/ranges/lookup on device in the production GPU path.
Small host metadata computation and unit extraction are acceptable; no host-built
per-atom cell lists or per-cell CPU scheduling loop disguised as GPU construction.

The default path updates geometry and rebuilds sorted cells each evaluation. It
does not need a materialized atom-pair list. Traverse surrounding cell ranges and
compute current distances directly. Surrounding cells need not form one contiguous
array: process their ranges without concatenating atom selections. Do not write
an N-by-M distance matrix, Boolean matrix or tuple-index matrix before evaluation.
Use squared distances to filter candidates before square roots where applicable.

Preallocate reusable record, sort scratch, range, parameter, force and reduction
buffers. Capacity growth is allowed when needed; distinguish steady-state reuse
from growth. No silent neighbor truncation. Host scalar reads/synchronization
must have an identified reason and be included in timing where applicable.

## 5. Geometry and ownership

Only single-image GPU geometry is supported. Initial scope: axis-aligned
orthorhombic periodic cells, fully open cells and mixed periodic/open directions.
Periodic lengths must be finite/positive and strictly exceed twice the search
radius. Search radius is cutoff for fresh cells and cutoff+skin for buffered lists.
Open generators may be zero/infinite as supported by the existing CPU boundary;
do not apply finite periodic-cell requirements to vacuum axes. Validate finite
positions and positive finite cutoff. Use minimum-image displacements only on
periodic axes. Recompute displacements from current geometry every evaluation.

Reject multi-image geometry when GPU is explicitly requested. If an automatic
backend policy is introduced, document CPU fallback and avoid silent fallback
inside explicitly timed GPU execution. Existing CPU multi-image support remains.
Skew cells may be explicitly unsupported initially. Do not implement them using
componentwise fractional rounding and claim general closest-image correctness.

Wrapped cell aliases must be deduplicated in this single-image mode: small counts
of cells per periodic direction must not cause the same candidate to be traversed
more than once. Exclude zero-shift self pairs by atom identity, not distance>0:
distinct coincident atoms require the potential's defined origin behavior.

Baseline uses full directed traversal. For each center, accumulate its full force.
Each physical pair contributes half its returned energy and virial in each
direction. The potential callback still returns unweighted physical pair values.
Apply the existing sign/index convention: d=Xj-Xi+H*s, g=(dE/dr)*d/r,
Fi+=g, W[a,b]-=d[a]*g[b]. Preserve the existing virial convention, not a newly
chosen transpose/sign. Validate whether a callback supports symmetric directed
evaluation; reject incompatible asymmetric/legacy semantics rather than assume.

## 6. Generic radial mathematics

Reuse the distance/species/output contract of NI.evaluate_pair! and its prepared
potential/scratch methods. GPU scalar views (one pair per lane/thread) can satisfy
the protocol without embedding CPU SIMD.jl vectors. Mathematical callback returns
E and/or dE/dr through requested output accessors; NeighbourList owns orientation,
cutoff masks, reductions, weighting and scheduling. No per-potential GPU neighbor
algorithms. PairPotentials owns GPU-compatible model preparation and tables.

Do not assume mutable CPU PairOutput or arbitrary CPU scratch compiles on device.
Implement suitable accessor-compatible local/output-buffer views and test actual
GPU compilation. Resolve preparation/device traits once, not by dictionary or
dynamic dispatch in hot loops. Stateless callbacks and explicitly supported local
scratch are sufficient initially; reject unsupported scratch/layouts clearly.

Support bundled radial models, including Lennard-Jones and SmoothPair where
applicable after inspection, and a custom device-compilable radial example. Reuse
their cutoff/shift/switch mathematics and derivatives. General legacy callbacks
need not all become GPU-compatible; leave CPU compatibility intact and report GPU
limitations. Anisotropic Cartesian callbacks are not required in this radial phase.
At r=0 require regular dE/dr=0 or explicit rejection; a safe denominator alone is
not permission to hide a singular potential. Inactive/padded geometry must be safe
before callback operations and inactive outputs must never be consumed.

## 7. Two execution strategies

### Baseline: central ownership
Assign each center a unique owner (one thread or cooperating subgroup/workgroup).
Traverse candidate cell ranges and accumulate only that center's force. Several
threads may reduce within a workgroup. No inter-workgroup force atomics are needed
when exactly one workgroup owns each center. Energy/virial use reusable partials
and parallel reductions. Benchmark thread-per-center versus cooperative execution
only where useful; do not multiply implementations without evidence.

### Experimental: candidate-thread tile
A workgroup owns a bounded central-atom tile. Each thread handles a candidate
atom from surrounding cell ranges and loops over central atoms, staging central
positions in shared/local memory for reuse. Split large cells/candidate regions
into bounded chunks; never require one cell to fit one workgroup. Candidate loads
should use contiguous sorted runs. Barriers must be reached uniformly, including
inactive tail lanes. Do not assume a fixed NVIDIA warp or AMD wavefront width.

Initially use the same full-directed/central-force-only mathematics as baseline.
Threads reduce separately for every central atom. If multiple workgroups split
one center's candidate region, use bounded tile partials and a separate reduction
(or an explicitly measured accumulation alternative). Unique force ownership
cannot be claimed simply because contributions are called central. Avoid a full
N-atom force buffer per workgroup. One workgroup iterating candidate chunks can
avoid global partials; choose and document the tested tradeoff. Half-list/Newton
two-endpoint evaluation is a later experiment, not required for phase completion.

Both strategies must use the same potential mathematics, precision and geometry
rules. Expose a strategy switch without changing physical potential definitions.

## 8. Optional skin reuse (implement last)

Default skin disabled: rebuild fresh sorted cells every call. Optional skin>0
builds conservative candidate atom/tile connectivity at cutoff+skin, saving
reference positions/cell and the sorting permutation. Between rebuilds update
positions in that fixed sorted order and follow stored connectivity, not stale
geometric cell adjacency. Filter all current interactions at the physical cutoff.
For tiles use conservative atom/bounding-volume inclusion, never just center
distance without a bound. No cached distance/displacement reuse.

Rebuild before evaluation when 2*maximum displacement>=skin, using a conservative
rounding policy. Cell/periodicity, count/order/identity/species, cutoff/skin or
precision changes invalidate preparation. Any cell change rebuilds initially.
Use documented unwrapped coordinates/tracked winding, or a rigorous conservative
fallback when motion cannot be established from input; do not silently infer
arbitrary unwrapped motion from wrapped positions. Refresh image geometry under
the documented input convention. Device reduction plus a small host rebuild flag
read is acceptable. No automatic skin tuning. If it does not improve speed, report
that result; it is not a prerequisite for the first working GPU milestones.

## 9. Explicit exclusions and delivery

No GPU ACE/MACE, GPU AD, MPI, electrostatics beyond existing short-range models,
matrix-operation distance algorithm or upstream AtomsCalculators proposal. GPU
AD requests must reject clearly rather than silently returning detached outputs;
plain inference need not raise merely because an AD package is loaded. Preserve
all existing functional CPU derivatives and do not advertise new differentiation.

Deliver source, tests, runnable examples, benchmark scripts, phase status and
reproduction commands. Describe support per tested backend, not generic claims
of CUDA/Metal/oneAPI support based solely on KA compilation intent.
