# CPU ACE specification

## Purpose and scope

Implement an independent evaluator for the supplied ACEpotentials model, validate
it against the same mathematical model, and investigate complete CPU force-stack
performance. Reuse existing neighbour/pair infrastructure without assuming that
pair ownership implements site-based ACE. Preserve APIs where feasible.

Required this phase: model reconstruction/import; serial and SIMD CPU ACE;
preallocated workspaces; full species-grouped environments; energy/forces/virial;
linear ACEfit fitting; Lux and ChainRules support; expanded pair derivatives;
angular/fusion performance experiments. GPU execution and MACE layers are later.

## Model, parameters, state and scratch

Use an explicit numerical core plus convenient calculator wrapper. Illustrative
roles, not unverified API signatures:
- Model: species labels, basis specifications, dimensions and structural indices.
- `ps`: authoritative trainable values (readout and exposed radial/mixing values).
- `st`: fixed numerical tables and nontrainable state.
- Environment: adjacency, image shifts, species buckets, masks, packing mappings.
- Workspace: mutable scratch, capacity and exclusive execution ownership.
- Calculator: owns model/ps/st, environment/workspaces and inference caches for MD.

Treat structural arrays as device-transferable later even if not trainable.
A parameter may be held fixed during linear ACEfit fitting while remaining in ps
for nonlinear training. Fixed coupling coefficients are not automatically trainable.
Use nested parameter structures/array shapes compatible with Lux; the existing
flat scalar/vector NamedTuple pair interface is not a universal restriction.
Keep fixed reference terms represented explicitly.

There must be one authoritative parameter set. Precontracted coefficients and
other derived numerical data are caches. Parameter updates refresh those caches
explicitly. Training differentiates through their construction or supplies the
correct adjoint to the original ps; never detach an old cache from new parameters.
Do not put reusable overwritten scratch in Lux state as though it were persistent
forward data. Use the actual Lux/LuxCore interfaces available in supported versions.

## Numerical geometry boundary

AtomsBase extraction and unit conversion occur outside hot loops. Core positions
are explicit plain numerical arrays in original atom order, with cell H using
columns. Compute d_ij = X_j + H*S_ij - X_i from current positions and fixed metadata;
keep spatial/packing permutations discrete but propagate numerical gradients
through gather/scatter. A species change or topology change invalidates metadata.

Existing API units remain eV, Å, eV/Å; public outputs retain documented units.
Species identity in the neighbour system is strict ChemicalSpecies identity.
Reference models may use element-only channels: define and test that conversion,
including a clear unsupported isotope/name policy rather than silently conflating
species. Positions and parameters must not be rounded to a fixed Float type in
an AD path. Fast ordinary Float32/64 paths can remain specialized.

This phase differentiates positions and parameters at fixed cell/topology.
Keep cell explicit for virial calculation and future extensions, but cell AD and
pressure-coupled trajectories are not required. A rebuild is a discrete operation;
ignore derivatives of selection/sorting/image choices, not derivatives of geometry.
Correct cutoff smoothness and valid topology are preconditions for derivatives.

## Complete environments

Build all directed neighbours for each centre. Do not apply half-list ownership
inside ACE site evaluation. Each centre's site energy is counted once. Distinct
periodic images with the same atom index remain distinct edges; exclude only the
zero-shift self-pair. Nonzero self-images must behave correctly for forces/virial.

Baseline: use existing snapshot image semantics to validate a full environment.
Production: reusable flat adjacency/offset arrays, exact shifts, mapping back to
sorted/original atom indices, species bucket offsets, lane masks and current
geometry. Avoid calling copying `neighbours` queries in the hot loop. Counting/
bucketing is sufficient; a general sort per site is unnecessary.

Pack accepted neighbours by species for uniform parameter selection. Record tail
masks; replace inactive inputs with values safe for each transform/normalization
before evaluation and mask values and tangents. Zero times NaN is not protection.
Use one conservative cutoff for preparation with per-component physical masks
when needed. Cutoff taper belongs to the imported radial/model definition once.
Cutoffs stay fixed parameters for now. Rebuild every call initially; persistent
capacities do not imply Verlet reuse. A later skin requires explicit invalidation
and separate tests before any reuse claim.

## Basis computation

Radial: reproduce imported transforms, envelopes, recurrence coefficients and
channel mixing. Implement scalar reference and SIMD across neighbour distances.
Joint values/d/dr recurrence evaluation should reuse intermediate arithmetic.
Keep only necessary recurrence orders live; store/consume results as produced.
Provide value-only and derivative request paths. Do not assume every model uses
a simple l-independent radial family.

Angular: initial backend is SpheriCart with exact reference harmonic conventions.
Pinned reference has SpheriCart 0.2.4; current main kernels inspected earlier may
be different. Verify actual selected APIs, gradients and normalization locally.
Allow normalization/sign/order conversions explicitly. Real/complex conventions
must be matched throughout coupling maps, not only at final harmonic output.
Provide a stable internal angular interface for later SIMD block backend.

Density projection: sum species-resolved radial-angular products for a complete
centre, respecting imported ordering. Materialized neighbour-major radial/angular
arrays are the transparent baseline. Explicitly test empty species groups.

Products and couplings: begin with selected product index lists, standard sparse
coupling maps and in-place dense multiplies as appropriate. No TensorOperations
core dependency. Preserve structural nonzeros independently of learned values.
Reverse products without division by inputs, so zeros are valid. Shared-product
DAGs and compiled nonzero execution lists are later measured optimizations.

Linear inference: where B=C*P and E=c'*B, precompute C'*c for a fast readout.
Retain original B/parameter coordinates for fitting and diagnostics. This algebra
only applies across fixed linear maps; do not move parameters through nonlinear
transforms or silently drop reference/pair terms.

## Forces, virial and composition

Use explicit position derivatives/reverse numerical kernels for ordinary force
evaluation. No finite-difference production forces. Edge adjoints g=dE/dd give
+g to centre and -g to neighbour with the displacement convention above.
For ACE, sum over directed site contributions without an arbitrary factor 1/2.
Pair terms retain their existing unordered-pair ownership.

Keep repository convention W=-sum(d ⊗ g); force is on the central endpoint.
Verify the existing transpose/sign convention using homogeneous strain finite
differences. Virial is configurational, in eV, with no volume normalization.
Stress conversion is a separate explicit convention.

Support sums of pair, ACE and fixed reference components. Shared environment
preparation is desirable when cutoffs permit, not required before correctness.
Honor AtomsCalculators additive ! output conventions and original atom order.

## Storage and concurrency

Preallocate model-dependent buffers and capacity-manage neighbour/batch-dependent
buffers. Grow outside inner kernels when needed. Reset all active regions; avoid
stale values after a smaller system/species group follows a larger one. Separate
borrowed workspace outputs from owning public results with explicit lifetime rules.

CPU scratch is per exclusive worker slot/task, not a global mutable buffer or a
thread-ID assumption across migration/yields. Keep one shared sorted force array.
For parallel ACE, compute full environment write footprints; a lower-cell-owned
pair schedule cannot be reused unchanged. Reuse colouring algorithms with new
footprints; tasks in a phase are disjoint and each task stays internally serial
initially. No hidden full N-atoms force array per worker.

Keep extraction/preparation allocations separate from warmed numerical-kernel
allocations and public result ownership. Bumper can remain an optional local
strategy for irregular temporaries. It is not the main ownership mechanism.

## Differentiation contract

Required for pair and ACE:
- energy -> positions and explicit parameters;
- forces -> positions and explicit parameters;
- virial -> explicit parameters.
Optional: virial -> positions; report unimplemented derivatives explicitly if
requested, rather than returning a false zero. Zero virial cotangents must not
prevent otherwise supported energy/force differentiation of a combined result.

Force reverse mode needs energy-Hessian actions and mixed position/parameter
products. Implement VJPs/JVPs without materializing full Hessians/Jacobians.
An rrule that returns forces for energy is insufficient by itself. Support the
explicit force-output primitive with rules/fallbacks and test it via the chosen
AD backend. A correctness-oriented scalar/local ForwardDiff fallback is acceptable
initially; label cost/scaling and avoid replaying whole-system evaluation once
per parameter as the final ACE training strategy.

Do not differentiate sorting, species, masks or graph indices. Pack/gather
numerical values and scatter-add their adjoints through retained maps. Treat
integer image choice as fixed locally. Repeated destinations, including self
images, must sum correctly in forward/reverse accumulation.

A pullback must own saved values or recompute from preserved forward inputs.
Neither parameters/positions nor workspace scratch may be read after caller
mutation unless the API explicitly requires immutability and enforces its lifetime.
Preserve the existing safe saved-pullback behaviour. Include repeated pullback
calls, reused workspaces and parameter update tests.

Lux should expose explicit ps/st and a practical value interface plus property
calculation interfaces. Keep in-place kernels behind derivative-aware boundaries.
Use Zygote as the first reverse-mode integration target (already a test dependency),
plus ForwardDiff/local analytic references where applicable. Do not imply support
for every AD backend or arbitrary differentiation order.

Validate a short deterministic fixed-cell MD rollout through forces, including
initial position/velocity and parameter sensitivities as appropriate. This is an
integration check, not an MD integrator package deliverable or long-time gradient
accuracy claim. Pressure coupling and virial position AD remain deferred.

## Fitting

ACEfit performs linear regression only with internal basis parameters held fixed.
Assemble total-energy features, negative position derivative force rows, and
virial rows under the established convention. Preserve basis/parameter ordering.
Subtract fixed baseline contributions from labels. Support weights and a documented
smoothness regularizer derived from basis metadata. Transform priors consistently
when changing basis coordinates; sparsity/normalization changes do not preserve
regularization automatically. Keep parameter vectors dense unless true structural
sparsity justifies another representation.

Place ACEfit integration in an optional extension or isolated training module.
Verify ACEfit's actual matrix-solver/assembly APIs. Dense CPU design matrix is
sufficient initially. Assemble and fit linear pair coefficients jointly with ACE
where supported. Nonlinear radial/pair parameters use Lux with a separate optimizer.

## Future GPU boundaries

Model data, explicit parameters/state, discrete graph metadata and owned scratch
must remain separable. Backend layouts and batch sizes may differ. Future device
workspaces will own persistent arrays; shared memory/register planning, asynchronous
lifetimes and force reduction require separate implementation. Do not equate CPU
SIMD lane width or CPU colouring with GPU execution. No GPU implementation now.
