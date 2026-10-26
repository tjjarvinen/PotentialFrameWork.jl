# Validation and acceptance

Use existing test/reference infrastructure. Add tests of mathematical and ownership
contracts; avoid tests merely asserting file/module layout or mirroring loops.
Run meaningful focused suites during development, then the existing full suite.
Use Julia 1.10-compatible and current supported environments when available;
never rewrite pinned reference manifests to force dependency resolution.

## Numerical behavior

- Pair and ACE energy/force/virial versus current independent serial/reference
  implementations, in Float32 and Float64 with appropriate existing tolerances.
- Forces versus finite differences of energy and virial versus homogeneous strain,
  including nonsymmetric geometries to detect transposed tensors.
- Pair ownership, energy normalization and legacy force-sign adapters.
- Zero atoms where existing APIs support them, isolated atoms, empty neighbor
  environments and coincident physical atom rejection where required.
- Orthorhombic/skew cells, open/slab/wire periodicity, unwrapped inputs, small
  periodic cells, multiple image edges and nonzero self images.
- A self-image can cancel atom force while retaining energy/virial: check this.
- Existing raw-array cell orientation and original atom ordering.
- Safe masked lanes and species padding; poison unused outputs where useful to
  reveal stale consumption. No NaN leakage from inactive lanes.

## Requests, units and parameters

- Every nonempty subset of Energy/Forces/Virial; reordered/duplicate public
  requests normalize consistently; invalid property requests fail clearly.
- Legacy EnergyOnly and EnergyForcesVirial and old callback wrappers still work.
- Lengths supplied in convertible units agree with Å input.
- Native callbacks see declared numeric conventions; explicitly unsupported
  model-unit conventions fail in preparation rather than producing wrong values.
- Unitful energy, forces and virial retain existing result types/conventions.
- Updating explicit parameters changes output without rebuilding model structure;
  no stale precontracted weights or bundled-parameter precedence bugs.
- Custom radial channels, source-backed invokelatest safeguards, masks, pruning,
  restricted bases and species-order projection retain existing behavior.

## Scheduling and batching

- batch_size 1, 2 and a larger nondivisor of the site count; invalid zero/negative
  values reject. Confirm direct single-site versus specialized batch dispatch.
- Cross batch sizes with representative SIMD widths (including width 1) and
  serial/threaded execution.
- Unequal neighborhood lengths/species groups and a smaller final batch.
- Selected centers nothing, empty, singleton and reordered selections; forces
  still reach nonselected neighbors. Preserve duplicate/invalid selection rules.
- Validate complete endpoint footprints, phase barriers and exclusive scratch
  ownership; no scratch indexing by Threads.threadid().
- Changed topology/capacity rejects stale prepared execution or is explicitly
  re-prepared at the public boundary.
- Returned public results survive later evaluation, resizing and workspace reuse.

## Interoperability

Use an independent simple backend with different storage/accessors and the same
potential callbacks. Include an analytic angular site energy with independently
checked derivatives, not just a sum of radial pairs.
Document its limited geometry capabilities and test explicit rejection outside
them. Native performance paths must not need this reference representation.

Retain existing NeighbourLists optional query tests and dependency-free behavior.
This phase does not require a production upstream neighbor-provider adapter.

## AD and external boundaries

Re-run existing ForwardDiff, ChainRules and Zygote coverage:
- Energy -> positions and parameters.
- Forces -> positions and parameters, including force-loss pullbacks.
- Virial -> parameters; preserve explicit unsupported derivative errors.
- Custom radial parameter derivatives and selected centers/projections.
- Delayed/repeated pullbacks after caller position/parameter mutation and
  workspace reuse preserve forward-evaluation inputs.
- Lux explicit parameters/unchanged state and AtomsCalculators calculate/get/set
  parameters, additive force destinations and requested-property output shapes.
- Fitting/design-matrix and model persistence regressions.

Do not silently reroute differentiable calls through Float32/64-only prepared
inference, mutate saved AD inputs, or claim new optimized AD.

## Performance and evidence

Before/after, measure warmed prepared pair and ACE evaluation separately from
public extraction/topology construction, with identical inputs/threads/options.
Inspect hot-loop inference and allocations using normal Julia tooling.
No per-site general allocation or abstract callback dispatch in native hot loops;
bounded prepared scratch/output capacity growth is allowed.
Retain the existing one-shared-force-array production policy.
Record batch-fallback overhead and slot-local gradient-buffer memory.
Batch fallback is a correctness feature, not evidence of throughput improvement.

Existing small reference structures suffice for regression checks; do not create
a major benchmark campaign or infer large-system speed from tiny timings.
Report noisy timings honestly and investigate clear regressions before completion.
Record exact versions/hardware/commands and any unavailable runtime checks.
Do not invent passing tests, unsupported GPU claims or speedups.

