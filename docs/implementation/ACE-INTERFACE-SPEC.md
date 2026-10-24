# ACE interface specification

## 1. Architectural and compatibility contract

Separate mathematical basis description and provenance, selection/transformation,
parameters, backend execution plans, environment and mutable scratch. Do not turn
plotting metadata into data traversed by hot kernels. Keep the root module thin;
place functionality within ACE or a narrowly scoped supporting submodule.

Keep species and parameters as runtime data and hot loops concretely typed.
Generalize concrete ResolvedModel/RadialChannel restrictions only where needed;
preserve existing public constructors/loaders through adapters. No required
ACEpotentials, plotting, CUDA or AMDGPU dependency in the core. Preserve practical
Julia 1.10.10 compatibility and the separate pinned reference environment.

Maintain original-order output, Å/eV numerical conventions, Unitful public
boundaries, full periodic images, additive AtomsCalculators output conventions,
one shared production force array and exclusive scratch ownership. Never retain
precontracted weights after authoritative parameters or structural plans change.

Every new selection/transformation must work consistently in scalar reference,
CPU execution, diagnostics, linear fitting, persistence and supported AD paths.
Avoid an inference-only feature that silently trains a different model.

## 2. Descriptor inspection and visualization

Expose owning descriptor metadata, summaries and arrays usable without plotting.
Record original identity/index, current identity/index, central readout channel
where relevant, radial/harmonic/product dependencies, invariant coupling terms,
correlation order(s), normalization conventions and coefficient correspondence.
Namespace provenance by source basis identity; array position alone is not a
stable identifier. Do not store a single tuple of angular momenta as if it fully
described an invariant that is a sum of coupled products.

Expose exact metadata where available. Mixed-order rows need their actual order
set/range. Nominal maximum body order is correlation order plus the center, not
a guarantee of distinct atoms. Current algebraic sum(n+l) is a defined index
weight; do not call it a physical polynomial degree for an arbitrary callable.
Unknown/undefined degree should be represented explicitly, not invented.

Add UnicodePlots as a weak dependency through a package extension, using its
direct API rather than requiring Plots.jl. Verify actual API and compatible
versions. Initial visuals: counts grouped by correlation order/index weight,
coefficient magnitude distribution, selected radial curves, coupling sparsity,
and a grouped count heatmap where informative. Plot functions return plot objects;
do not automatically plot from show or model loading. Support bounded selections
and readable labels. Coefficient magnitudes are exploratory, normalization- and
representation-dependent, not an automatic importance criterion.

Advanced design-matrix correlation/SVD analysis and new per-descriptor energy
contribution APIs are deferred. Existing site_diagnostics remains available;
add no new public site_energies function or AtomsCalculators property.

## 3. Basis masking, pruning and restricted construction

Use one immutable, inspectable selection description for complete invariant
features: IDs/indices and available order/index-weight criteria. Function names
in earlier discussion are illustrative; choose a small coherent public API.

Masking keeps original basis dimensions, ordering and reversible selection data.
Enforce disabled coefficients as effective zeros on every evaluation and after
updates, with zero sensitivities for disabled trainable entries. Masking must
remain enforced during fitting/readout installation, not be a one-time assignment
an optimizer can undo. Preserve the unmasked parameters for reversible changes.

Structural pruning removes unwanted complete invariant features and unused
dependencies. Trace coupling rows -> structurally needed product columns ->
density projections -> radial/harmonic requirements. Keep shared dependencies.
Do not prune based on learned coefficient cancellations or temporary zeros.
Compact indices and transfer numerical tables/parameters consistently; include
nonlinear radial scales, pair scales and fixed reference terms. A compact angular
index map may coexist with a backend computing full harmonics through lmax;
report such remaining work honestly. Return explicit coefficient/dependency maps
and provenance. Fully empty retained invariant sets must be supported.

Restricted construction applies this selection while constructing a numerical
model/execution plan from an existing basis specification/import. It must not
construct a full evaluator/workspace first merely to discard it. It may read the
existing full specification. A general new ACE invariant/coupling basis generator
is NOT part of this phase. Document this boundary clearly.

Unselected pair/reference components are preserved by ordinary invariant pruning.
Expose explicit component inclusion policies when requested; never silently drop
the baseline. Same retained coefficients: masking and pruning give equivalent
E/F/W, not necessarily the predictions of the original untrimmed model.

## 4. Central-index selection

Support fixed central ranges/vectors at explicit numerical and calculator entry
points, and carry the same choice into fitting/AD. Default: every center. Indices
refer to original system order, independent of internal permutations. Preserve
input selection order in metadata; summation order need not imply output layout.
Reject duplicates and invalid/out-of-bounds indices. Empty selection is valid and
returns zero E/F/W with correctly shaped full-system forces.

The energy is sum(E_i for i in selected centers), including enabled reference,
pair and ACE site components. Preserve the current imported directed pair-site
convention; selecting centers does not select complete unordered pairs. Each
selected center uses all relevant neighbors. Accumulate derivatives onto ALL
affected atoms; do not slice the system or zero forces outside selected centers.
Filter/rebuild CPU execution phases using full force-write footprints. Full-system
environment preparation is acceptable initially; report its costs separately.

No dedicated central-species basis selector now; users can form indices by species.
Fixed center selection must remain explicit discrete metadata during AD. Changing
regions, transition forces and adaptive-resolution physics are deferred.

## 5. Neighbor-species correlation-order selection

Accept one selected neighbor species or a species group and an inclusive integer
range L:K of selected-species factor counts. A group has ONE combined counter,
not independent counters. Central identity never contributes to this count.
Validate species identity using the model's existing element-channel policy.
Reject malformed ranges/unknown species; document empty-group behavior. Order
above a descriptor's degree contributes zero; no-selected-neighbor cases remain
well-defined. Unrestricted selection must recover the original model.

For each density index a split A_a = U_a + V_a, where V sums neighbors belonging
to the selected species group and U sums the rest, with identical radial functions,
scales, harmonics and cutoffs. For each product, compute polynomial coefficients
of product_a(U_a + t V_a). Keep coefficients k in L:K, then apply the unchanged
angular coupling map and readout. This is a transformed feature, not a claim that
existing mixed-species descriptor rows carry definite chemical labels.

Use truncated coefficient convolution; no symbolic expansion into all species
assignments. Retain intermediate orders needed to reach the maximum, even when
the minimum is high. Reverse without division by U/V/product values; zeros must
be safe. Reuse radial/angular evaluation for U and V; preallocate polynomial and
reverse scratch in CPU plans. Benchmark costs without claiming selection is
necessarily faster. Preserve complete invariants, rotation and permutation symmetry.

Factor count is correlation order, NOT number of distinct neighbors; repeated
atom terms are included. Default: sum retained k components with the ORIGINAL
shared coefficient. Do not silently introduce independent trainable coefficients.
Preserve k provenance so a future explicit independent-coefficient construction
can extend this; independent-component fitting is not required in this phase.

Reference energies and nonselected-neighbor pair terms have selected order 0.
Selected-neighbor pair terms have selected order 1. Include/exclude them according
to L:K as well as explicit component enable flags. Thus L>=1 cancels references;
0:0 retains nonselected terms. This is filtering the algebraic decomposition of
the chosen site energy, not an isolated subsystem or a unique physical partition.

Composition of center selection, invariant masking/pruning and species-order
projection must be supported. Preserve original coefficient mappings and make
operation semantics unambiguous. Avoid a separate exponential-size species basis.

## 6. Custom radial callable interface

Support functions and callable structs, with explicit configuration and trainable
parameters separate from model structure and scratch. The core protocol should
write preallocated basis values and optionally joint distance derivatives; a
convenience adapter may accept allocating functions. Document costs of adapters.

Metadata must declare output dimensions, ordering, radial/angular correspondence,
central/neighbor channel dependence, cutoff, envelope ownership and any algebraic
index weights used by priors. Do not infer dimension by calling at an unsafe
distance. Allow n,l-dependent functions rather than assuming all l share a basis.
Select species channels outside inner loops where possible. Preserve imported
Agnesi/spline functions through the same protocol and their existing fast path.

Values and derivatives must use promoted generic numeric types. Scalars and
SIMD.Vec inputs represent independent distances with the same output channel
ordering. Avoid Float64 conversions, concretely typed scalar-only outputs and
host allocations in the SIMD fast protocol. Provide an analytic joint evaluator
and a documented scalar AD fallback where supported; do not promise automatic
SIMD differentiation of arbitrary functions. It is acceptable to adapt scalar
AD derivatives lane-wise as a labelled fallback. Keep nonlinear parameters in ps,
including custom parameters, and differentiate them through the evaluation.

For force-loss training, analytic distance derivatives alone are insufficient:
their position/parameter derivatives must be available through generic AD-compatible
code or explicit derivative products. Keep existing E->X/ps, F->X/ps and W->ps
contracts; preserve existing optional W->X support where it already works. Report
unsupported capabilities clearly before returning incorrect gradients.

Write a portable example with nontrivial trainable parameters, smooth cutoff,
scalar/SIMD implementations or one generic implementation, and independent checks.
Adapt an ACEpotentials radial object in the optional pinned reference environment
as an interoperability example, verifying actual API, normalization, cutoff and
ordering. Radial compatibility alone is not full model compatibility.

## 7. Persistence

Extend the existing save_calculator/load path; do not invent a conflicting second
authoritative parameter tree. Preserve loading old native/resolved formats. Use
explicit schema versions and migration/validation for additions. Save model,
selection policies, coefficient/dependency maps, custom radial definitions and
authoritative parameters/state. Rebuild environment/schedule/scratch/caches.

Support three custom radial routes:
1. Registered built-in/extension implementation: versioned ID and data config.
2. Explicit Julia source: source text/file(s), entry point/factory, interface
   version, configuration, dependency/version metadata, source hash and numerical
   reference evaluations. No attempted recovery of source from arbitrary closures.
3. Opt-in tabulation into a supported self-contained numerical spline representation.

Custom-source loading defaults to disabled and requires explicit allow_code=true
or an equivalently clear API. Check this BEFORE evaluating ANY saved source.
Evaluate definitions in a fresh namespace; it avoids name collisions, not malicious
code. Handle Julia world-age correctly at the loading boundary, without putting
invokelatest in the inner radial loop. Do not silently install packages or execute
saved shell commands. Report missing dependencies/unsupported versions. Reference
checks are validation, not a security sandbox; a source hash establishes identity,
not trust. Self-contained source must not depend on captured session globals.

Validation covers parse/load/entrypoint, applicable methods, dimensions, values,
derivatives and declared numeric capabilities. Reference checks should be
type-tolerant with documented tolerances. Source-only scalar success is not proof
of SIMD, AD or GPU support. Saving an unsupported callable without a source or
registration should fail clearly and leave existing files intact.

Tabulation is explicit approximation, never a silent replacement. Document domain,
interpolation, boundary/cutoff behavior and value/derivative validation error.
Sample away from training nodes too. Do not assume arbitrary callables have the
existing Agnesi transform; define a suitable tabulation mapping. Tables are a
fixed snapshot and do not preserve derivatives with respect to original custom
parameters. For force training/vibrations report second-derivative smoothness
limitations; do not claim equivalence beyond checked tolerances.

## 8. Fitting, regularization and documentation

Design-matrix columns use selected/transformed features in their authoritative
coefficient coordinates. Fixed baselines obey component and species-order rules.
Mask/prune mappings must also govern fit installation and nonlinear parameters.
Preserve Unitful label handling. With c_old = T*c_new, transfer a prior operator
as R_new = R_old*T; do not assume it stays diagonal for general maps. Undefined
custom radial index weights need an explicit user prior or a documented neutral
choice; never silently assign misleading smoothness degrees.

Documentation must distinguish CPU scalar, CPU SIMD, training AD and FUTURE GPU
capabilities. GPU requires device-compilable operations, concrete device-compatible
state, no host-only calls or unsupported allocations; scalar-per-thread radial
evaluation is a natural GPU path. Accepting SIMD.Vec does not establish GPU support.
Explain masked branches/safe inactive lanes, origin/cutoff policy, units, output
ownership and derivatives. GPU requirements are documented, not runtime-tested
support in this CPU-only phase.

Defer changing regions, central-or-neighbor rules, independent group counters,
general basis generation, composite-model compilation, QM/MM/ML/MM physics,
new site-energy API, long-range electrostatics, GPU/MACE and major benchmarks.
