# Development milestones and acceptance gates

Implement M0–M8 in order, with small coherent changes and progress recorded in
STATUS.md. Independent contract/import tasks can overlap conceptually; there is
no instruction to delegate to subagents. Fix failures before optimizing the next
stage. Gates are verification checkpoints, not mandatory permission pauses.

## M0 — Reinspect and reproduce reference

Read the current local source/test tree and data/ace_ref. Record actual versions,
fixture hashes, component conventions and baseline test results. Preserve core
Julia compatibility; use the supplied reference environment separately. Load the
saved model, never refit by accident. Re-evaluate on reloaded XYZ coordinates and
create full-precision additional references where useful. Validate force provenance.

Gate: original evaluator reproduces a documented rounded/full-precision reference
relationship; the existing suite passes, or concrete pre-existing failures are
identified without suppressing them. Update current status above historical notes.

## M1 — Explicit numerical inputs and pair derivative proof

Introduce the smallest shared interfaces needed for explicit X/ps/st/environment
and scratch ownership. Preserve pair callback and calculator APIs through adapters.
Keep the existing SIMD inference path working. Prototype position and parameter
AD on SmoothPair and a nonlinear pair example. Extend explicit parameters beyond
energy-only requests. Implement force-position and force-parameter VJPs and virial-
parameter differentiation, initially with a verified correctness fallback if needed.

Do not attempt AD through `prepare!`'s wrapped Float arrays or assume the existing
energy-only rrule supplies higher derivatives. Define original/sorted mappings,
nonzero-image handling and pullback lifetimes now. Nested parameter shapes can be
introduced without rewriting all models into Lux layers at once.

Gate: finite-difference derivative-product checks, actual Zygote integration and
saved-pullback/workspace-mutation tests pass. Existing callers remain valid.

## M2 — Resolve and import the saved model

In the pinned reference environment inspect the constructed model and export its
complete numerical basis specification and parameter mapping into a versioned
resolved representation. Implement a narrow parser/adapter for this representation
and an import convenience path from the provided recipe JSON. Capture implicit
pair basis, Eref and actual angular/coupling conventions. Use explicit assertions
for component dimensions and supported cases. Keep ACEpotentials out of repeated
numerical evaluation. Record licenses/attribution for any reused code.

Gate: raw-to-resolved parameter mapping and basis metadata match reconstruction;
radial/angular samples agree; unsupported schemas/components fail clearly.

## M3 — Full environment preparation

Build serial complete directed environments with exact periodic image labels and
species groups. Use the current snapshot semantics/oracle for reference and flat
capacity-managed adjacency for the production route. Pack neighbours into uniform-
species SIMD blocks; retain all endpoint/image maps and masks. Compute numerical
displacements from current X. Rebuild each public call initially.

Gate: full directed image multiset matches the independent oracle on small/skew/
multi-image/open cells; padding and changed inputs are safe; every production
force destination can be enumerated. Test nonzero self-image derivatives explicitly.

## M4 — Serial independent ACE evaluator

Add the ACE functional submodule. Implement scalar radial reference, SpheriCart
adapter, species density projection, selected products, sparse invariant maps,
readout and all reference/pair contributions. Add explicit forces and virial using
reverse numerical stages. Keep intermediate diagnostics available outside hot paths.
Integrate AtomsCalculators energy/force/virial requests and owning outputs.

Gate: intermediate samples and final predictions agree with the original evaluator
on all 33 fixtures using geometry/precision-appropriate tolerances. Independent
energy-position and strain finite differences, invariance and zero-coefficient
cases pass. No call to ACEpotentials evaluates production outputs.

## M5 — Linear fitting and composition

Provide structure-level basis values and force/virial design rows. Integrate the
actual ACEfit API in an extension, CPU dense assembly first. Implement label
selection/weights, fixed reference subtraction, basis smoothness metadata and
coefficient updates. Support fixed and jointly fitted linear pair contributions.
Keep inference precontraction and original training coordinates equivalent.

Gate: design matrix times coefficients equals independently evaluated E/F/W;
synthetic coefficient/prediction recovery tests pass; small energy+force fitting
and virial assembly work; updated ps refreshes inference caches correctly.
Rank deficiency need not permit unique coefficient recovery: test predictions and
identifiable coefficients appropriately.

## M6 — SIMD radial, workspaces and CPU parallel execution

Use SIMD radial recurrence across neighbours; joint values/derivatives; species
parameters bound outside lane arithmetic. Keep a scalar baseline. Capacity-manage
numerical buffers and environments, separating public preparation/output costs.
Add fixed-map/readout precontraction where valid, with tested cache invalidation.

Build new full-environment footprints and phase schedule; reuse colouring logic,
not the existing pair footprints. Allocate exclusive scratch slots independent of
thread IDs. Accumulate into one shared sorted force array. Profile scheduling and
species packing overhead before adding more elaborate strategies.

Gate: SIMD/precision/thread agreement, actual footprint tests, reuse/growth/reset
checks and allocation reports pass. Record generated SIMD arithmetic where relevant.
Report performance including complete preparation; no speedup requirement.

## M7 — Complete Lux/ChainRules integration

Finish shared Lux adapters and pair/ACE rules. Required outputs and tangent inputs
are listed in ACE-SPEC. Implement efficient reverse numerical stages and local
higher-derivative products needed by force pullbacks; use recomputation/saved data
with explicit ownership. Support differentiable exposed radial/mixing parameters,
not only linear readout. Verify precontraction gradients into original ps.

Gate: real Zygote energy/force/virial-parameter losses, position force VJPs, local
ForwardDiff/finite-difference references, optional-virial handling and short fixed-
cell rollout sensitivities pass. No full Hessian is required. Clearly document
any remaining optional unsupported derivative and backend limitations.

## M8 — Angular block and fusion experiments

Retain SpheriCart. Implement/test an explicit SIMD block angular evaluator using
compatible recurrence/conventions where practical; verify actual versioned source
and license. Adapt generic constants, normalization, derivatives and result types
rather than assuming SVector{3,Vec} is a supported public API.

Compare standalone SpheriCart batched, generated single-input (if available in
selected version), SIMD angular blocks, and fused radial-angular-density kernels.
Tile to limit live intermediates; do not assume all harmonics/gradients/projection
accumulators stay in registers. Consider recomputation for force contraction.
Benchmark multiple angular cutoffs, neighbour/species counts, SIMD widths and
values/gradient workloads. Distinguish neighbour lanes from optional centre lanes
for later contractions. No automatic backend selection without repeatable evidence.

Gate: numerical agreement and complete E/F/W performance results, including packing
and allocations. Keep reference backend selectable; slower experiments can remain
optional, be removed with evidence, or be documented as unsupported if a verified
technical blocker exists. Do not fabricate a successful fused implementation.

## Separate later phase — GPU and MACE

After the CPU phase, plan device data, kernels, batching, persistent allocations,
force reduction and AD. This handoff does not authorize implementation of GPU/MACE.
No acceptance gate here depends on CUDA/AMDGPU hardware or dependencies.

## Completion report

Update README/examples, STATUS.md and validation/benchmark documentation. Give
commands, versions and limitations. Preserve first-phase benchmarks as historical.
State which backends and derivatives actually work; report all incomplete required
gates explicitly. Avoid a universal speedup claim from one TiAl fixture. Do not
silently replace unimplemented functionality with a reference-package call.
