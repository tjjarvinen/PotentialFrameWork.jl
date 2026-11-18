# Generic site pullbacks and equivariant ACE training — Codex instructions

## Objective and authorized scope

Implement this work in the user's local VS Code checkout of [PotentialFrameWork.jl](https://github.com/tjjarvinen/PotentialFrameWork.jl). Run Julia and validation locally. This document was prepared by source inspection; it does not claim new test results.

Deliver two sequential, independently reviewable parts:

| Part | Deliverable |
| --- | --- |
| A | A model-independent `site_vjp!` protocol and serial chunk backward executor, demonstrated with an independent nonlinear model, scalar ACE E/F/V derivatives, and first-order equivariant readout training. |
| B | Complete equivariant ACE force-loss training through a nonlinear invariant energy head, including mixed derivative products, head-parameter derivatives, and end-to-end validation. |

Complete Part A's acceptance gate before proceeding to Part B. Both parts are authorized; do not stop for another approval after Part A. Keep their changes and evidence clearly separable. Resolve routine signature and implementation choices from the actual checkout and record the resulting contracts.

The neighbor-list implementation and interface must serve models beyond ACE/MACE. Keep ACE mathematics and irrep conventions outside generic execution code. Keep NeighbourInterface and AtomsCalculators as protocols; execution, allocation, reductions and scheduling belong to implementation modules in PFW. Do not extract or publish new standalone packages in this work.

Out of scope: persistent Channel/task execution, NUMA placement, Dagger, GPU AD, GPU equivariant kernels, MACE message passing and derivatives into incoming atom features, trainable topology/species/cell/cutoffs, model-format redesign, new fitting-data requirements, and a broad performance optimization campaign. Preserve existing pair/GPU behavior and checkpoint/import compatibility. A small invariant head is a validation/example component, not a MACE implementation.

## Inspect the checkout first

Planning baseline: commit `6d91cf5f5c6d7514816040c40a0bbc2286626737`, “Neighborlist interface improvement (#18)”, 2026-11-17. Build on subsequent local changes; do not reset to this commit. Inspect applicable AGENTS.md files, branch, status and uncommitted changes. Preserve user files, reference models, manifests and historical measurements. Do not delegate to agents unless the user or applicable repository instructions explicitly authorize it.

Read the current versions of:

- `docs/site_chunks.md`, `docs/ace_features.md`, `docs/ace_interface.md`.
- `docs/implementation/STATUS.md`, `CPU-Site-Chunk-Codex-Instructions.md`, `ACE-Equivariant-Codex-Instructions.md`, and relevant ACE/neighbour interface specifications and validation documents.
- `src/submodules/NeighbourInterface/NeighbourInterface.jl`, `src/submodules/NeighbourList/chunks.jl` and local site/batch protocol implementations.
- `src/submodules/ACE/{chunks,serial,features,calculators,radials}.jl` and relevant product/angular code.
- ChainRulesCore, ForwardDiff, SpheriCart and LuxCore extensions.
- Existing site-chunk, scalar AD, feature AD, custom-radial, interface and rollout tests.

Existing functionality to reuse: functional scalar `ace_evaluate` has E/F/V pullbacks into explicit parameters and positions. `ace_vjp_impl` uses a displacement directional differentiation of analytic reverse stages. `ace_features` and `feature_vjp` support first-order positions/parameters for arbitrary output blocks, but do not establish general force-loss training through composed nonlinear heads. Chunk evaluation currently supplies inference and indexing-transpose tests, not a complete AD rule. Verify these facts against the working tree before refactoring.

## 1. Mathematical and ordering contract

For one central site i, use displacements `r_ij = X_j - X_i + H*s_ij`. A local model may return site energy E_i, edge displacement gradients g_ij, and optional per-site features h_i. For conservative scalar energy, `g_ij = ∂E_i/∂r_ij`.

The executor's current assembly convention is:

```text
E += E_i
F_i += g_ij
F_j -= g_ij
W -= r_ij * transpose(g_ij)
```

Use the Frobenius inner product for matrix cotangents. Given global E/F/W cotangents, the executor supplies:

```text
energy_bar_i = energy_bar
g_bar_ij = force_bar_i - force_bar_j - transpose(virial_bar) * r_ij
direct_r_bar_ij = -virial_bar * g_ij
```

The model VJP returns the displacement derivative of its local output contraction. Add `direct_r_bar_ij` when positions are requested, then scatter `-r_bar_ij` to i and `+r_bar_ij` to j. Energy/gradient/output contributions must be combined without duplication. Handle absent requests and zero cotangents explicitly. If another existing backend uses a different virial tensor orientation, adapt at its boundary; do not copy its transpose convention blindly.

Gather per-site feature cotangents by central atom and preserve model-defined block layouts. Convert original-order output cotangents into cell order exactly once at an appropriate boundary. Preserve individual periodic-image edges while accumulating aliased atom destinations. Self-image endpoint position contributions cancel; their feature/parameter/virial contributions need not vanish. Inactive padding contributes nothing and must never access invalid destinations.

Topology, cell and species are fixed for differentiation, and topology must remain valid near the differentiated geometry. Use current numeric positions plus saved image shifts when recomputing geometry. Do not differentiate frozen forward displacements in place of the position argument.

## 2. Part A — generic local derivative protocol

Introduce `site_vjp!` as a declared model callback in the appropriate pure protocol module. The illustrative signature is:

```julia
site_vjp!(input_bars, parameter_bar, prepared_model,
          site, parameters, output_bars, scratch; request)
```

Adapt argument order to repository conventions. Before using it widely, document its exact signature, public buffer types, axes, accessors, ownership and accumulation semantics. Do not leave the interface defined only by examples or by ACE-specific implementation details.

### Required local semantics

- Compute the VJP of the whole site map, not independent edge calculations. Neighbor coupling within a site is allowed.
- Local output cotangents cover energy, edge displacement-gradient outputs and model-defined per-site feature outputs, independently or combined. Use explicit requested-output contracts rather than assuming all models provide every field.
- Initially, differentiable inputs are displacements and explicit parameters. Incoming atom-feature derivatives are reserved for MACE; no implementation is required now.
- Support parameter-only, displacement-only and combined input derivative requests. Do not compute a position Hessian product merely because a parameter-only force pullback was requested.
- Accept arbitrary numeric feature cotangents; no ACE L/parity or magnetic-axis knowledge belongs in NeighbourInterface/NeighbourList.
- Define which callback destinations accumulate and which are overwritten. Prefer additive cotangent destinations: executor resets fresh pass/slot accumulators; callbacks add only active contributions. Document return value and alias restrictions.
- Borrowed site/output/scratch views cannot escape a callback. A model may recompute local intermediates or use an explicitly owned saved context, but must not retain mutable worker scratch as a backward tape.
- Normalize ChainRules zero/thunk cotangents at the AD boundary or through a documented adapter, so every model need not duplicate their handling. Preserve the distinction between a mathematical zero and an unsupported derivative.

### Capabilities and parameter accumulation

Expose a small explicit way to query/validate supported local output-to-input derivative requests. Distinguish an energy/feature VJP from a VJP of gradient outputs. A supported first derivative must not imply force-loss or arbitrary-order AD support. Unsupported active requests fail clearly before partial results are published; zero cotangents must not force irrelevant unsupported work.

Provide allocation/reset/add operations for model-compatible parameter accumulators. Preserve nested NamedTuple field names, real scalar leaves and array shapes. Permit model specialization without requiring flattening as the public protocol. Any internal flattening must have a stable validated inverse and must not be repeated per edge merely for reduction convenience. Immutable scalar leaves may require returning/rebinding an updated tree or mutable internal accumulator leaves; specify this precisely.

Keep these operations serial now but suitable for later exclusive worker accumulators and a final reduction. Do not allocate one full parameter gradient per site or store per-edge parameter Jacobians. Include empty parameter trees.

Define backward batch compatibility. A documented fallback that invokes `site_vjp!` on active sites of the existing batch is sufficient. Optimized SIMD/batched derivative kernels are not required. Do not imply that batched forward inference automatically establishes batched AD.

## 3. Part A — serial backward executor and adapters

Add a numerical differentiable entry point with explicit positions and parameters. Supply a custom ChainRules rule around this mathematical operation; internally run chunk forward/backward passes. Do not attempt to make arbitrary `foreach_worker!` closures automatically differentiable. Execution state, topology, derived caches and model tables are not trainable parameters.

The backward executor owns output-cotangent permutations, center gathers, E/F/W-to-local cotangent conversion, endpoint scatter, parameter reductions and original-order export. Use the existing bounded chunk destination mapping. Preserve selected-center semantics: a selected site can contribute to unselected neighbor atoms. Keep per-site feature buffers in cell order between passes and export only at the public boundary.

Expose the supported AD path naturally through the relevant numeric/Lux and explicit-parameter calculator entry points. Test the intended public paths, including units at the AtomsCalculators boundary. Do not merely add a callable local helper while leaving advertised training APIs on an undifferentiable mutable route. A functional reference remains available, but do not silently use it as the entire backward executor and claim chunk integration is complete.

Backward calls must not mutate caller parameters. Preserve safe delayed/repeated pullbacks after caller X/parameter/model/topology mutation, later forward calls and workspace reuse. An owning numerical snapshot followed by local recomputation is acceptable. Avoid blindly retaining/copying an entire executor with scratch, stale cache and historical fields; document saved data and costs. Distinct explicit states must remain independent; one state has one active invocation. Failures must not publish partial backward results.

### Independent model

Implement a small analytic nonlinear site model in test/example support, independently of ACE. Include a nonlinear parameter and preferably a simple cross-neighbor term so that treating edges independently would fail. Validate energy, gradient-output and optional simple feature VJPs against independent calculations. For example, a model with `q = sum(abs2, r_j)` and `E = a^2*q/2 + b*q^2/4` has nontrivial parameter/displacement derivatives and coupled neighbors. Ensure any generic adapter is not specialized on ACE types.

### Scalar ACE

Extract/reuse the existing local scalar reverse and mixed-product mathematics behind `site_vjp!`. Preserve radial/pair/readout parameter authority, selected-center behavior, supported projection policies and reference contributions. Avoid duplicate derivative implementations that can drift. Maintain the existing functional path as an independent execution reference. Reuse the displacement directional method where suitable; no full Hessian or system Jacobian should be formed.

### First-order equivariant ACE

Adapt existing feature forward/VJP logic to site/chunk execution, preserving arbitrary output irreps, parity, species-resolved mixing, radial scales and supported nested callable parameters. Existing output arrays have `(magnetic component, output channel, original atom)` axes; retain that public convention and output ownership. Generic feature storage must not assume a fixed SVector width or a single block.

Demonstrate direct training against scalar/vector/tensor feature targets using first-order derivatives. Keep the current invariant-head force evaluation working. This part does not yet claim force-loss training through that composed head.

## 4. Part A acceptance gate

Before Part B, pass locally:

1. Independent local VJP finite-difference/directional checks, including arbitrary feature cotangents and a cross-neighbor contribution.
2. End-to-end E/F/W local-to-global adjoint identities and gradient checks. Use nonsymmetric virial cotangents to expose transpose mistakes.
3. Scalar ACE chunk-versus-functional gradients for energy, forces, virial and combined losses, into every supported parameter subtree and positions.
4. Parameter-only, displacement-only and combined requests; verify via instrumentation where practical that omitted derivative work is not executed.
5. Equivariant first-order feature and parameter gradients, non-scalar blocks, single active blocks, parity choices, multispecies mixing, supported custom radials and a direct-feature training example.
6. Nontrivial atom permutation, forced chunk splitting, inactive SIMD lanes, batch size 1 and >1, selected centers, repeated images, self images, empty/isolated systems and zero density factors.
7. Missing/zero/thunk cotangents, empty parameter trees, scalar leaves, unsupported active capabilities, repeated and delayed pullbacks, caller mutation and state reuse/failure behavior.
8. Relevant existing regression suites with unchanged scientific tolerances. Preserve Float32/64 forward support; state precisely which backward precisions are supported. Feature execution need not acquire a new Float32 mode.

Record this gate and proceed to Part B without asking for permission.

## 5. Part B — higher derivatives for general equivariant training

For feature blocks `h = Φ(X, θ)` and a differentiable invariant head `E = ρ(h, α)`, support parameter gradients through `F = -∂E/∂X`, including both θ and α. Position derivatives through forces are also required for the established fixed-cell derivative contract; arbitrary third/higher derivatives are not implied.

The calculation must include both the changing feature Jacobian and the changing head feature cotangent. Holding `∂ρ/∂h` fixed during force-loss differentiation is incorrect for nonlinear heads.

Implement the feature JVP and mixed/Hessian products needed to compose these operations. A useful identity for a force cotangent B is:

```text
<B, F> = -D_X E[B]
D_X E[B] = <∂ρ/∂h, J_Φ,X B>
```

The derivative of this contraction must propagate through both factors and into head parameters. Choose explicit ChainRules rules at the required feature reverse boundary, an equivalent differentiable JVP construction, or another validated local product formulation. Document exactly which nesting/composition is supported. A first-order `rrule(ace_features, ...)` alone is insufficient, and mutation inside `feature_vjp` must not be left to unsupported nested reverse AD.

Reuse the existing scalar directional differentiation strategy where suitable. Support arbitrary numerical cotangents in every block and avoid constructing full feature-position/parameter Jacobians or full Hessians. Local parameter AD for custom radial callables remains acceptable; replaying the entire feature calculation for each parameter is not the intended algorithm.

Custom radials need mixed distance–parameter derivatives for force-parameter gradients and second distance derivatives for force-position gradients. Validate generic arithmetic, nested AD tags, analytic derivative consistency and supported values-only adapters. Fail clearly for declared unsupported callables. Fixed radial tables/coupling/cutoff metadata remain nondifferentiable.

Integrate these products with Part A's model-independent callback/executor. Do not add a second ACE-only system backward loop. Higher-order composition may need an additional JVP/local rule; it should complement `site_vjp!`, not overload a first-order capability flag to promise arbitrary nested AD.

## 6. Part B acceptance and training examples

Use small generated systems and models; no new fitted reference potential or user dataset is needed.

- Include multiple scalar, vector and tensor channels, supported even/odd parity, multiple species and nonzero active non-scalar blocks.
- Validate a smooth nonlinear invariant energy head. Norm-squared block contractions with trainable smooth nonlinear combinations are suitable. Include trainable head parameters with nonzero force sensitivity; an additive energy bias alone cannot test this.
- Check energy/feature first derivatives and force-loss parameter gradients independently for each mixing block, radial scales, supported nested custom radial parameters and head parameters. Use central finite differences over a sensible step range and directional checks, not only agreement between two AD paths sharing the same rule.
- Check feature JVP/VJP adjoint identities and selected mixed/Hessian products. Test force-position products against finite differences and symmetry identities where applicable.
- Test energy/force combined losses and optional direct feature-target losses. Check appropriate rotational/inversion transformation of outputs, forces and derivative contractions; physical losses must use an invariant head, while arbitrary output cotangents remain legal for VJP checks.
- Repeat essential periodic-image, lifetime, zero-block and chunk/batch checks from Part A on the composed model.
- Provide a deterministic, small executable training example through the supported public/Lux path. Use generated teacher targets, nontrivial initial parameters and an energy/force objective; optionally include direct equivariant targets in their declared frame. Report finite gradients and reproducible loss reduction with conservative optimizer settings. Do not claim coefficient recovery for a non-identifiable parameterization.
- Re-run the relevant regressions after these changes, including existing scalar force-loss and short rollout tests. Add a short composed-feature rollout check if needed to substantiate force-to-position differentiation.

Use existing testing conventions and justified numeric tolerances. Do not weaken existing gates to accommodate errors.

## 7. Storage, performance and compatibility

Keep worker scatter storage bounded by chunk capacity and numerical scratch bounded by maximum site/batch sizes. Final atom gradients, retained features and global parameter accumulators naturally scale with their respective outputs. Do not allocate a full atom-gradient array per site/worker or a parameter gradient for every retained site. Saved snapshots can scale with topology/inputs; report them separately from reusable scratch.

Include small warmed measurements separating preparation, forward, pullback construction/snapshot, backward execution and output export. Report memory and allocations for parameter-only versus combined derivatives and functional versus chunk paths. These are diagnostics, not a required speedup claim. No optimized backward SIMD kernel or broad benchmark campaign is required.

Preserve explicit parameter authority on every forward/backward call. Derived inference coefficients must not bypass differentiation into their source parameters or become stale. Keep current units, output ownership, model/checkpoint formats and supported CPU pair/GPU interfaces. Prefer existing optional AD dependencies and preserve current core Julia compatibility. Resolve test extras in separate local environments when needed; do not rewrite user/root/reference manifests merely to run tests.

## 8. Documentation and final handoff

Document the precise public protocol and provide a short recipe for a new model author: local outputs, derivative capabilities, parameter accumulator hooks, `site_vjp!`, optional JVP and registration with numerical execution. Explain that forward callbacks alone do not imply force-loss support. Include single-site and batch-fallback examples.

Update STATUS.md and relevant neighbour/feature/AD documentation with separate Part A and Part B records. Report:

- Final signatures, buffer/axis/accumulation contracts and changed public entry points.
- Supported output-to-input derivative combinations, precisions and AD composition orders.
- How functional, chunk, Lux and AtomsCalculators paths select execution and preserve saved inputs.
- Exact local validation commands, Julia/dependency versions, passed/failed gates and focused numerical evidence.
- Training configuration and loss before/after, plus allocation/saved-storage diagnostics.
- Any limitations or required follow-up, without claiming unvalidated GPU, threading, MACE or arbitrary-order AD support.

Finish with a concrete reviewable implementation and evidence for both parts. If a gate cannot be completed, identify the precise blocker and retain passing work; do not substitute a functional-only fallback or a first-order demonstration for the requested chunk force-loss/equivariant capability.
