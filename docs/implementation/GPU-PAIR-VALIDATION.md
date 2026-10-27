# GPU phase 1 validation and measurement

## Correctness references

Use the existing CPU evaluator and independent small O(N^2) single-image sums.
For Float32 compare against CPU evaluation on the same converted geometry and
parameters as well as the original higher-precision reference. This distinguishes
input rounding from execution/reduction differences. Do not demand bitwise identity
for floating-point reductions. Set and document justified absolute/relative tolerances
per model/type/scale; include max per-atom force error and normalized energy error.

Required checks:
- Float64 input -> Float32 compute; Float32 input -> Float64 compute; both native
  types. Assert prepared tables, intermediates and outputs have the selected type.
- Analytic force checks against energy differences and virial against strain
  differences at suitable step sizes, away from nonsmooth boundaries. Test cutoff
  boundary behavior separately. Check total force for isolated pair systems.
- Every supported Energy/Forces/Virial property subset, normalized/empty requests,
  sign conventions and full-directed energy/virial factor. No unrequested output use.
- Atom permutation invariance, correct original-order force restoration, same
  species aligned after sorting, species-dependent lookup, missing species errors,
  integer species round trips, and explicit parameter updates after preparation.
- Empty, one-atom, isolated-atom and distinct coincident-atom cases; singular-origin
  rejection versus supported regular origin; safe inactive/padded lanes.
- Open and orthorhombic periodic/mixed cells; positions outside the box by several
  periods, boundary crossings, few-cell wrapped aliases, empty cells and sparse
  clusters separated by vacuum. Check cutoff eligibility and rejection at equality.
- Explicit rejection of multi-image and unsupported skew GPU requests; CPU support
  must remain unchanged. Separate neighbor-image eligibility from physical cutoff.
- Arbitrary cell populations, central tile tails, candidate tails, multiple workgroups
  sharing a center, and strategy equivalence. Check race-sensitive cases repeatedly.
- Workspace growth, reuse, precision/config changes, stale topology prevention and
  result independence after later calls. No silent truncation or offset overflow.

## Skin trajectories

Test motions below and at/above the rebuild threshold; two atoms approaching from
outside cutoff but inside list radius; newly interacting pairs after rebuild; boundary
crossings and the documented winding convention; motion returning toward reference;
cell/species/order/cutoff changes; skin=0. Compare each step to fresh construction.
For tile candidates include a case where center-based pruning would miss a corner
atom interaction. Check floating-point margins conservatively. Test ambiguous motion
handling explicitly rather than claiming wrapped coordinates reveal winding.

## Regression and dependency checks

Run existing CPU pair, interface, periodic-image, species, calculator and AD tests.
Do not alter historical CPU numerical/ownership contracts to accommodate GPU storage.
Run CPU-only package loading without AMDGPU/CUDA installed/loaded. Opt-in real-GPU
tests must clearly distinguish unavailable/skipped from passing. In the dedicated
GPU validation command, an expected but unusable device must fail visibly.

Keep KA CPU-backend tests as useful kernel checks, but label them separately from
real-device results. No requirement to test a vendor without available hardware;
report it as untested. Verify the project's supported Julia version in dependency
resolution; do not confuse the separate reference model environment with core support.

## Benchmark structures

Generate with AtomsBuilder using verified installed APIs, then convert to an
AtomsSystems system using its supported constructor. Generate/replicate/convert/
perturb once before timing. Use fixed seeds and report exact atom count, density,
cell, periodicity, species fractions, cutoff, skin and perturbation amplitude.

Include approximately 10k and 100k atoms (replication need not yield round counts),
single-species argon LJ and multispecies cases with an appropriate bundled model.
Use physically reasonable spacings and bounded perturbations avoiding overlaps.
Include nonuniform occupancy or a slab/open case to expose load imbalance. Small
reference systems are for correctness, not evidence of GPU throughput.

## Timing protocol and reporting

Warm up compilation, preparation and each measured specialization. Use backend
synchronization at measurement boundaries; do not time only asynchronous enqueue.
Separate one-time allocation/compilation, capacity growth and warmed steady state.
Do not force synchronization between every kernel solely for routine timing.

Report two public workloads:
1. AtomsSystems -> owning unitful result, including extraction, transfer, sorting,
   search, interactions, reductions, output permutation and host result copies.
2. Prepared device-numerical execution, including its normal sorting/search and
   reductions, without convenience host output copies. State whether input upload
   is included and how positions are updated.

For skin mode separately report rebuild steps, reused-list steps, actual rebuild
frequency and amortized time over an evolving trajectory. Static repeated inputs
alone are not an adequate skin benchmark. For default mode include fresh-cell
construction every call. Compare strategies using matching property requests and
precision, including E-only and full E/F/W; compare CPU results on the same inputs.

Break down extraction/upload, key generation, sorting/unpacking, cell ranges/lookup,
candidate construction (when present), interaction evaluation, reduction and output.
Use selected profiling runs if per-stage synchronization distorts total execution.
Record buffer/device memory footprint and steady-state allocation behavior. No full
force array per workgroup or dense distance matrix should hide in scratch.

Record hardware/software versions, commands, repetitions and variation. Report
both wall time and useful work; candidate-pair rate alone can reward excess work.
Sorting may remain negligible, skin may not help, and candidate tiling may lose:
all are valid measured outcomes. No prescribed speed target is a correctness gate.
