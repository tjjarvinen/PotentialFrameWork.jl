# Implementation milestones

Proceed in order, using small reviewable changes. Gates are checks, not permission
pauses. Complete this phase before proposing GPU implementation. Refine internal
types as source inspection warrants, preserving the stated mathematical behavior.

## I0 — Baseline and API inventory

Record current branch/hash, Julia/dependency versions, available hardware and local
reference hashes. Read latest STATUS plus existing ACE tests and extension code.
Run the existing test suite in the supported environment. Preserve the pinned
ACEpotentials Manifest and never refit/overwrite the supplied model accidentally.
Identify method dispatch tied to ResolvedModel in calculator, AD, Lux, fitting and
persistence before introducing wrapper/generalized model types. Sketch concrete
types and a short public API in documentation; proceed without another design pause.

Gate: reproducible baseline or precise environmental/pre-existing blockers recorded.

## I1 — Basis metadata and complete-feature selection

Add backend-independent descriptor descriptions, basis summaries, stable provenance
and coefficient/dependency mappings. Define correlation-order versus nominal body
order, mixed-order rows, index weights and unknown custom degrees. Implement an
immutable selection specification, coefficient masking, structural pruning and
restricted construction from existing/imported specs. Integrate parameter updates
and preserve fixed pair/reference components unless explicitly selected otherwise.
Support empty invariant bases and repeated transformations without losing provenance.

Gate: pruning/masking equivalence on E/F/W and parameters, complete invariant
symmetry, no dependence on current numerical coefficient values, actual reduction
in appropriate structural arrays. No promise of proportional runtime reduction.

## I2 — Central indices

Thread fixed central indices through numerical evaluation, calculator calls,
diagnostics where applicable, design matrices, Lux and derivative rules. Implement
all/range/vector/empty behavior, validation and original-index semantics. Preserve
full neighborhoods and all force destinations. Filter phase tasks without violating
complete write footprints or allocating an N-atoms force array per worker.

Gate: selected sum and its derivatives agree with independent site sums, including
nonselected-neighbor forces and imported pair/reference accounting.

## I3 — Species correlation-order projection

Implement the scalar U+tV polynomial oracle, retained-order sum and analytic reverse
without division by factors. Keep coefficients shared. Add explicit transformed
metadata and component-order policy. Integrate selection with masks/pruning and
center subsets. Implement reusable CPU polynomial/reverse scratch and projection
buckets; reuse radial/angular evaluation. Extend fitting and AD coherently.

Gate: explicit small-factor expansion, sum-over-orders reconstruction, full-minus-
excluded presence identity and E/F/W/AD equivalence. Shared coefficient installation
and derivative mapping verified. No exponential species expansion in production.

## I4 — Radial protocol and custom callables

Introduce the minimal protocol and typed metadata/config needed for arbitrary
callables, with existing RadialChannel as a supported implementation. Preserve its
SIMD optimized path. Add values-only and joint derivative paths, explicit custom
ps, analytic derivatives and a documented AD fallback. Exercise scalar and SIMD
inputs, masks, Float32/64 and mixed parameter precision. Generalize affected model
and extension dispatch rather than bypassing AD or the numerical evaluator.

Add a smooth portable custom example and a pinned-environment ACEpotentials radial
adapter example. Verify the real reference API before relying on it. Custom
channel replacement must validate compatibility with existing radial/density maps;
it is not automatic generation of new invariant coupling tables.

Gate: unchanged imported predictions; custom E/F/W, parameter/position and force-loss
products correct; unsupported custom capabilities fail clearly. No GPU implementation.

## I5 — Persistence and approximation export

Extend checkpoint schema with old-format reading preserved. Include selections,
provenance/maps and custom basis config with one authoritative parameter tree.
Implement registered radial serialization, explicit source bundles and opt-in
tabulation. Do not serialize mutable runtime caches, compiled method objects or
arbitrary closures. Refuse code loading before evaluation unless opted in; isolate
namespaces and handle load-time world age. Never install dependencies automatically.

Gate: fresh-process/relocated round trips, code-disabled nonexecution, corrupted
metadata rejection, missing-dependency diagnostics and atomic save behavior.
Tabulation reports independently checked approximation errors and fixed-parameter
snapshot semantics. Preserve historical resolved and calculator checkpoint loading.

## I6 — Optional plotting and usage documentation

Add UnicodePlots weakdep/extension and small reusable plotting functions consuming
ordinary summary arrays. Document exact APIs with runnable examples: inspect,
mask/prune, restricted construction, centers, species order, custom radial basis,
save/load and plot. Provide a CPU/SIMD/AD/future-GPU requirements table plus a
validation helper returning separate capabilities and diagnostic reasons.

Gate: core loads without UnicodePlots; extension plots actual data, empty/selected
bases behave sensibly; no screenshot/golden rendering tests. No invented model data
presented as measured. GPU status remains unverified/deferred.

## I7 — Integration and completion

Run focused checks while developing, then full appropriate supported-version suites
and optional-extension checks. Record small measurements of kernel allocations,
structural sizes and species-order overhead, not a major benchmark campaign. Finish
README/examples and STATUS with reproduction commands, actual outcomes and limits.
Run git diff --check. List any incomplete gate explicitly rather than silently
substituting an unsupported implementation. Do not begin GPU work in this phase.
