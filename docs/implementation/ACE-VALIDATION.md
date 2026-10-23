# Validation and measurement requirements

## Baseline and environments

Run the existing suite before/after material changes. Preserve current supported
Julia compatibility; test a compatible 1.10 environment and a current supported
environment where dependencies allow. Reference reconstruction uses its supplied
Julia 1.13.2 Manifest separately. Do not overwrite that Manifest to force core
compatibility. Optional extensions must not break ordinary package loading when
their weak dependencies are absent.

Record commands, Julia/dependency versions, hardware, thread count and actual
results. Do not describe recorded historical test counts as freshly verified.
If runtime/network is unavailable, report which checks could not run.

## Reference and import

- Verify the exact saved model hash; never accidentally refit the oracle.
- Check all 270 flattened parameters' mapping into actual model components.
- Match species channel ordering, radial definitions and implicit defaults,
  normalization/sign/order, product/coupling indices, Eref and pair contributions.
- Compare radial/angular samples, density coefficients, products and invariants
  where accessible, then all 33 fixtures' E/F/W. For internal basis changes compare
  through explicit mappings rather than asserting identical arrays blindly.
- Recompute reference outputs on reloaded rounded XYZ geometries for tight tests.
  Store full precision or compare live. Use absolute and relative tolerances suited
  to Float32/64, output rounding and accumulation order. Large Eref values can hide
  errors in relative total-energy checks: also compare residual/component energy.
- Verify model predictions are distinct from original dataset labels. Corrected
  per-atom `force` fields are the prediction fixtures; energy and virial use eV,
  forces eV/Å. Add tests that prevent the old system-property export mistake.

## Numerical and geometry checks

Use independent original-order all-image oracles on small systems, not a second
wrapper around production geometry. Include:
- open/full/mixed periodicity, skew cells, tiny multi-image cells, wrapped and
  unwrapped positions, nonzero self-images and exact image reconstruction;
- duplicate atom indices with distinct image shifts, reverse edge directions,
  no zero-shift self edges, isolated/empty systems, empty species buckets;
- masks/tails, safe dummy directions before normalization, cutoff-edge values and
  documented origin policy. Physical coincident atoms differ from padded lanes;
- strict ChemicalSpecies mapping versus model element channels;
- force conservation, translation and permutation equivariance/invariance, proper
  rotations of energy/forces/virial where the model is invariant;
- scalar/SIMD widths 1/representative 2/4/8/16 and Float32/64 as supported;
- zero coefficients and zero density/product factors; no derivative division by zero;
- energy-only calls and additive AtomsCalculators ! outputs; output ownership;
- changed positions, cell, species, atom count, model/parameter shapes and cutoff
  policy; growing and shrinking workspace capacities and stale-tail clearing.

Compare forces to energy finite differences and virial to homogeneous strain of
positions AND cell. Preserve repository convention: for strain epsilon,
  dE/dt = -sum(W .* transpose(epsilon)).
Include shear as well as dilation, and component contributions to expose signs
or pair/site counting mistakes. Avoid numerical derivative stencils at nonsmooth
cutoff/image-selection boundaries unless testing their documented behaviour.

## AD products

Test each requirement for pair AND ACE, using actual backend calls and independent
finite differences on small well-conditioned fixtures. For ACE include several
parameter blocks, not only readout, if radial/mixing values are trainable.

1. grad_X E = -F and parameter energy gradients.
2. For a fixed force cotangent B, define L=sum(B .* F).
   Compare its reverse gradients into X and ps to finite differences of L.
   This exercises Hessian and mixed-derivative actions without full matrices.
3. For a fixed virial cotangent C, compare grad_ps sum(C .* W).
4. Combined E/F/W requests with zero virial cotangent must support the required
   E/F derivatives. Optional unsupported W->X with nonzero cotangent reports an
   explicit error/NotImplemented as appropriate; it must not be a silent zero.
5. Masked/inactive lanes contribute zero tangents without NaN/Inf intermediates;
   packed scatter gradients sum repeated indices correctly.
6. Compare native parameter path, precontracted inference path and their adjoints.
7. Saved pullbacks survive workspace reuse and caller input/parameter mutation
   according to the owning contract; repeated calls are consistent. Include
   zero/thunk cotangents and tangent projection/shape checks.
8. Compare serial/threaded derivatives where provided. Threaded forward correctness
   is not evidence of a race-free reverse scatter.

ForwardDiff energy Hessians do not establish higher-order reverse-mode support.
An energy rrule returning -F does not establish force-output differentiation.
Document the verified backend and derivative order precisely.

## Fixed-cell MD integration

Use a small deterministic short velocity-Verlet or equivalent test rollout. Keep
velocities/masses explicit and avoid stochastic thermostats. Compare derivatives
of a terminal position/velocity objective with respect to initial positions and
model parameters; include initial velocity sensitivity if exposed. Choose a smooth
fixture, then test safe list rebuild changes separately where feasible. This is a
potential/AD integration test; no new production integrator is required. Do not
infer long-time stability or NPT differentiability from this test.

## Fitting

- A*c reproduces independent total energy, force and virial predictions including
  correct component subtraction/addition and original atom ordering.
- Synthetic linear fits recover identifiable parameters or predictions in rank-
  deficient cases; do not assert unique coefficients without full rank.
- Test energy+force fits, virial rows, nonuniform weights, missing labels and fixed
  baselines. Verify prior/regularizer scaling in the actual solve API.
- Linear pair and ACE feature columns compose correctly.
- Readout installation updates all inference caches; nonlinear parameter changes
  invalidate the right basis-dependent data and do not reuse stale design matrices.

## Parallel writes and memory

Enumerate actual ACE task destinations on small cases and verify they are subsets
of advertised complete-environment footprints; phases have disjoint write cells.
Exercise tiny periodic grids, skew multi-images, repeated destinations and uneven
populations. Repeat numerical comparisons at 1/4 and other useful thread counts
using separate Julia processes. Do not require bitwise equality without a specified
reduction order. Exclusive scratch must not rely on task thread-ID stability.

Measure force-array bytes versus N and workers. Production force storage must not
scale as N times worker count. Distinguish owning outputs, temporary edge adjoints,
saved AD data and model features from force-array copies. Report legitimate edge/
feature memory rather than hiding it under the one-force-array claim.

Warmed prepared numerical kernels should reuse allocated capacity; measure rather
than assume zero allocations. Preparation, public result conversion, thread task
creation and pullback saves may allocate and must be reported separately. Do not
make brittle timing assertions in CI. Use bounded allocation checks only where a
stable ownership guarantee justifies them.

## Performance experiments

Reuse the existing benchmark style but add an ACE-specific suite. Measure:
- extraction, topology/schedule, full environment/species packing;
- radial, angular, projection, products/couplings, readout;
- reverse/force/virial accumulation, output permutation;
- complete public E/F/W and energy-only calls;
- prepared unchanged-geometry diagnostics, labelled as such;
- allocations, model/workspace memory, lane occupancy, phase occupancy/imbalance;
- AD forward/backward time and saved-state memory separately where measured.

Vary system size, species mix, density/neighbour count, lmax/basis size, precision,
SIMD width and threads. Keep compiled warmup outside measured samples, avoid
concurrent benchmark jobs, record BLAS threads/Julia target/affinity, and report
variation. Control nested SpheriCart/BLAS threading when comparing outer parallelism.
Compare against the pinned ACEpotentials evaluator on the same model and geometry;
record whether setup/list rebuild and units/output conversion are included in both.

Angular experiments need identical standalone outputs and complete force timings.
A fused density kernel may win end-to-end while losing a standalone microbenchmark.
Inspect generated code/register spill evidence where practical; source-level named
locals do not prove register residency. Retain explicit backend choice initially.
No fixed speedup requirement, automatic crossover heuristic or GPU claim follows
from these experiments.
