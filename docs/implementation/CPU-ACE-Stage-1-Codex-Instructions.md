# CPU ACE optimization — Stage 1 Codex instructions

Date: 10 November 2026
Repository: https://github.com/tjjarvinen/PotentialFrameWork.jl

## Mission and authorization

Implement the first CPU performance improvement stage in PotentialFrameWork.jl. Preserve the existing implementation as a reproducible performance baseline and correctness reference. Establish larger, matched benchmarks against the established modern ACEpotentials evaluator, then implement an unrestricted-product fast path. Preserve general equivariant features, parity, explicit parameters, Lux integration, and all currently supported derivatives.

The longer-term research goal is to test whether PFW can outperform PACE while retaining broader flexibility. Beating PACE is not a stage-one acceptance requirement. PACE is optional and must not block the mandatory PFW/ACEpotentials work.

This document authorizes CPU changes and benchmark dependencies needed for this stage. It extends repository instructions only where needed for this scope. Preserve unrelated behavior and GPU functionality. Do not spawn subagents unless the user separately authorizes delegation. Proceed autonomously through routine reversible choices; do not stop at a plan.

## 1. Inspect the actual starting tree

Read AGENTS.md, the latest STATUS.md entries, CPU ACE/equivariant/interface specifications, and the GPU milestone-1 implementation record. Inspect actual code before deciding signatures or dispatch. Do not overwrite local models, manifests, uncommitted changes, or historical measurements.

Relevant starting points include:

- docs/gpu_ace_cpu_simd_audit.md
- benchmark/gpu_ace/audit.jl and audit_workload.jl
- benchmark/gpu_ace/results/2026-11-09-julia1.13.3-5eef1c1-dirty/
- benchmark/gpu_ace/export_legacy.jl, oracle_legacy.jl, oracle_tial.jl
- src/submodules/ACE/cpu.jl, projection.jl, interface.jl, features.jl
- src/submodules/NeighbourList/environments.jl and execution.jl

The checkout inspected during planning was ea076d3c18d5ae66776f1fa785e729f2eadfba54. This is context, not an instruction to reset the user's checkout. Pin the actual clean baseline or precisely record dirty source hashes and preserve a reproducible snapshot before editing.

The earlier SIMD audit used 128-atom serial prepared workloads. Its conclusions about coloring/thread concurrency do not generalize to larger systems. Initial radial assembly dumps inspected a default-argument wrapper; use audit-corrected artifacts for actual radial/SpheriCart evidence. Static instruction counts are not runtime percentages. Public Ti/Ti–Al–V analytic radials run lane-wise in the existing CPU callable adapter.

## 2. Preserve a usable baseline

Before optimization, record source identity, dependencies, model hashes, geometry hashes, and reproduction commands. Retain a separate baseline checkout/environment or equivalent immutable source snapshot. Do not duplicate the whole repository into the source tree.

Keep the existing product algorithm runnable as a reference path. Add a narrowly scoped explicit algorithm choice, such as reference/optimized/auto, using names consistent with the existing API. Tests must exercise each path. Default unrestricted inference may select the optimized path after validation; constrained cases must retain the reference/general algorithm. Preserve existing constructor and call compatibility.

A product-algorithm switch in the modified checkout is useful for controlled experiments, but it does not replace the immutable historical baseline: shared infrastructure changes can otherwise move both results.

## 3. Mandatory benchmark suite

Use BenchmarkTools and reuse existing acquisition, provenance, workload, and report machinery where practical. Add a CPU-specific entry point rather than requiring GPU packages or devices. Keep incompatible external environments isolated.

### Models

1. Retain the small TiAl reference model and its modern ACEpotentials reference.
2. Generate bounded medium and large models in the established modern ACEpotentials implementation and import identical models into PFW. Cover single-species and multispecies cases, increasing angular/correlation/degree complexity. Inspect current APIs; do not assume a constructor name identifies its backend. In particular, modern ace1_model need not execute old ACE1.
3. Fitting is not required. Deterministic, well-scaled nonzero coefficients are acceptable for performance fixtures; label these as synthetic and avoid physical-accuracy claims. Reject fixtures with nonfinite or severely ill-scaled predictions.
4. Retain the published Ti and Ti–Al–V models as supplementary legacy cases. Label ACEpotentials 0.6.6/ACE1 references explicitly, separately from modern ACEpotentials.

Record actual radial/angular/density/product counts, product-order distribution, coupling dimensions/nonzeros, species, cutoffs, and radial representation. Nominal invariant counts are not sufficient to compare computational work. Use structurally meaningful coefficient/species distributions; document active coefficients per central species where feasible.

If one modern fixture cannot be imported without major new representation work, finish compatible fixtures and explain the precise limitation. Do not substitute a different potential and call it matched. Keep generation reproducible; cache models with hashes, attribution, and licensing metadata where needed.

### Workloads and presets

- Approximately 128, 1,024, and 10,000+ atoms; deterministic structures, species assignments, and modest perturbations, with valid periodic cells for the cutoff.
- Float32 and Float64 in PFW, and matching precision in external evaluators where actually supported. Unsupported precision is explicit, not silently promoted.
- Serial measurements first, followed by representative process/thread scaling up to available physical cores. Fix and report BLAS threads and affinity where practical. Do not multiply every parameter combination unnecessarily.
- Quick preset: bounded smoke/correctness coverage plus representative small and 10k+ timing cases. Full preset: broader model/size/precision/thread coverage with configurable limits.
- Run expensive cases sequentially. Time budgets and sample counts must be configurable. Aim for at least 10 warmed samples for representative reported comparisons; retain actual counts, and flag under-sampled cases. Avoid impractically large default runs.

Report three distinct boundaries where available:

1. Prepared numerical combined E/F/V, excluding model/geometry preparation and owning output conversion.
2. Complete public calculator calls including their actual preparation and output costs.
3. Repeated-frame execution with an explicitly documented existing connectivity/rebuild policy.

Do not introduce CPU skin reuse in stage one merely to improve trajectory timing. If connectivity reuse is unsupported, report fresh rebuilds. External timing boundaries that differ must be labelled rather than presented as pure kernel speedups. Separate compilation/startup from warmed execution.

### Equivariant coverage

Add representative feature workloads with multiple scalar, vector, and tensor channels and supported parity choices. Include feature evaluation and selected existing parameter/position pullbacks. These are preservation and performance baselines, not a claim that the scalar optimization accelerates them. External feature comparisons are optional and require genuinely equivalent outputs/conventions.

## 4. Profile and instrument before selecting later work

Add opt-in diagnostic timing or isolated warmed microbenchmarks for:

- extraction and unit conversion;
- neighbor construction, species packing, footprint/color scheduling, and executor preparation;
- geometry/radius handling, radial evaluation, angular evaluation;
- forward density contraction;
- product values and reverse density adjoints;
- energy/readout accumulation, force back-contraction, force/virial reduction;
- owning output construction.

Keep production hot loops free of timers by default. Stage measurements can perturb cache behavior and must not be summed as proof of undisturbed total performance. Pair them with whole-call results and sampling profiles when available. Retain raw observations and allocations. Report scratch payload separately from allocation peaks.

Inspect generated code for concrete relevant Float32/Float64 signatures, including the actual callee and explicit parameters. Check unexpected widening, scalar dependency chains, bounds checks, indexed access, and allocation. Do not assume Float32 must provide 2x whole-evaluator speedup. Width 4 in both precisions is not a precision-optimal comparison by itself; retain configurable width experiments without implementing a new batching architecture.

## 5. First optimization: unrestricted product evaluation

Study product_projection! and its callers to preserve their exact output/adjoint conventions. The current implementation supports polynomial species-order projection through U and V and uses forward/reverse polynomial tables for each product.

When policy.orders === nothing, the unrestricted expression for a product with indices a_t is:

    P_j = product_t (U[a_t] + V[a_t])

Use a division-free prefix/reverse-sweep algorithm, or an equivalent stable linear-work algorithm per product, to compute values and weighted adjoints. Do not divide total products by factors: zero factors and repeated indices must work correctly. Reuse bounded scratch; avoid per-product allocation and unnecessary full polynomial-table clearing.

Species selection partitions U and V; unrestricted evaluation must remain correct with nonempty selections, not only the default V=0 case. Propagate derivatives consistently to both arrays as required by the existing interface. Preserve additive accumulation for repeated indices and multiple products. Seed reverse contributions using current weights exactly once, and avoid stale gradients between sites/calls. Keep masked/readout parameters authoritative on every invocation.

Constrained species-order policies continue through the existing general implementation. Do not infer mathematical equivalence for a restricted policy merely because selected values happen to be zero. Do not broaden dispatch based on untested conditions.

This optimization is initially for the prepared CPU scalar E/F/V path. Preserve the functional AD and general feature implementations; they may remain separate. Do not silently route AD through mutation or detach parameters. If optimized primitives are shared with those paths, validate every affected derivative contract explicitly.

Only perform additional small local changes when profiling supports them and their effect can be isolated. Do not roll graph construction, broad contraction rewrites, new reductions, or new batching into this first patch.

## 6. Optional PACE and spline comparison

The user confirmed that a PACE-enabled static LAMMPS installation exists in their development container at /opt/lammps-static/. Default candidate executable:

    /opt/lammps-static/bin/lmp

Accept a CLI/configured executable override and optionally an environment variable. PATH is not required. Probe the actual executable and inspect available packages/styles. User confirmation does not prove that this execution environment can access that path.

Use the isolated legacy ACEpotentials 0.6.6/ACE1 export route for the published Ti/Ti–Al–V models. It exports radial splines for PACE; confirm the actual exporter and interpolation semantics before relying on any historical node count. Preserve all pair, reference, unit, species-order, and energy-offset contributions. Include any separately required pair table in the LAMMPS input.

Maintain two clearly labelled comparisons:

- Analytic original: ACE1 versus PFW with exact analytic radial recurrences.
- Splined export: PACE versus PFW using equivalent exported spline data where supported. Include ACE1's splined evaluator only if verified usable.

Separate spline approximation error against the analytic original from implementation disagreement on the splined model. Matching nodes alone does not establish identical interpolation, derivatives, endpoint, or cutoff behavior. Verify whether export/lookup changes interpolation again.

A narrow optional spline importer is allowed if it can be implemented and validated without major redesign. Otherwise complete modern/legacy analytic CPU comparisons and record the PACE limitation. No full legacy-to-modern invariant-basis converter, refit, ETACE migration, or large dependency build is required for stage one.

Before timing PACE, verify energy, forces, and virial/stress conventions on multiple configurations. Test product and recursive modes when supported. Start with one process/core; report actual parallel backend rather than assuming OMP_NUM_THREADS parallelizes recursive PACE. Keep a LAMMPS process alive across repeated evaluations and exclude process startup. Distinguish fixed-geometry evaluation timing from MD trajectory timing, including neighbor rebuild and integration work. Never compare PFW prepared kernels directly to total container/process wall time.

If unavailable, incompatible, or numerically unverified, emit a structured skipped/unsupported result with the reason. A launched PACE run that fails correctness is a failed comparison, not a passing skip. Continue independent mandatory work and report all failures. Do not fabricate PACE timings.

## 7. Correctness and regression gates

### Product-level tests

Compare reference/optimized values and weighted U/V adjoints for both precisions, unrestricted policies with empty/nonempty selected-species sets, supported product orders, repeated indices, one/multiple zero factors, zero weights, mixed signs, and repeated calls with changed inputs/weights. Handle degree-zero/empty cases only if the actual model contract permits them. Use representative well-scaled and challenging finite inputs. Add independent finite-difference checks on a small subset, not only reference-algorithm comparisons.

### Full model tests

- E/F/V agreement with original PFW and independent ACEpotentials references.
- Force/strain finite differences on selected small non-symmetric structures; unit and virial-sign checks.
- Restricted-order projection, masking/pruning, central selection, custom radials, units, supported geometries, and parameter updates remain correct.
- Equivariant outputs preserve rotations, translations, parity/inversion, channel conventions, and existing feature/readout AD.
- Existing energy-to-position/parameter, force-to-position/parameter, and virial-to-parameter contracts continue passing. Preserve existing virial-to-position support or explicit unsupported behavior.
- Existing GPU and CPU-only loading behavior is retained. Run available relevant GPU regression checks if shared code changes; otherwise explicitly state unavailable hardware coverage.
- Preserve one shared production CPU force array, scratch ownership, safe colored scheduling, output lifetime, and thread safety.

Retain precision-appropriate established tolerances; choose and justify new tolerances before interpreting performance. Do not loosen existing contracts to obtain passing benchmarks. Distinguish Float32 implementation error against a rounded reference from application-level sensitivity against Float64.

Inference optimization does not imply faster AD. Document which path each derivative uses and ensure no stale state after parameter changes. Run meaningful existing package tests and optional integrations in their compatible environments; preserve supported Julia compatibility rather than silently raising it.

## 8. Reporting and acceptance

Write machine-readable results, a concise report, and exact reproduction commands. Record baseline/optimized source identities, Julia/dependency/compiler versions, hardware, external executable identity and version, thread/process/BLAS counts, model/geometry hashes, algorithm/backend, radial representation, timed boundary, samples/spread, host allocations, scratch payload, errors, and skip/failure reasons.

Show incremental speedup from the unrestricted-product change and total speedup against the preserved historical baseline. Report losses as well as wins. Never mix ACE1 and modern ACEpotentials labels or compare different potentials as matched. Preserve historical benchmark artifacts unchanged. Update docs/implementation/STATUS.md with real completed gates and remaining limitations.

Stage-one completion requires:

1. Reproducible baseline and quick/full CPU benchmark entry points.
2. Matched modern ACEpotentials comparisons including larger structures and meaningful model-complexity coverage, with precise documentation of any genuinely blocked fixture.
3. Stage diagnostics identifying remaining costs.
4. Correct unrestricted-product fast path, selectable reference path, and validated constrained fallback.
5. Preserved general equivariant and derivative contracts.
6. Actual performance evidence. If the fast path does not improve representative product-heavy workloads, investigate and explain; do not promote a regression or manufacture a speedup gate. Beating PACE or achieving a fixed Float32 ratio is not required.
7. Optional PACE results or explicit nonblocking omissions.

## 9. Deferred work and handoff

Do not implement shared-product DAGs, SIMD across central sites, broad density/force-contraction redesign, new CPU reduction ownership, CPU skin reuse, ETACE migration, MACE layers, or new GPU execution in this stage. Preserve boundaries for them.

Finish with a measured recommendation for the next stage: product sharing/species-specific plans if product work remains dominant; contraction/layout changes if those dominate; preparation/scheduling if full-call overhead dominates. PACE-inspired graphs and SIMD across same-species sites remain leading hypotheses, not predetermined winners.

Do not copy current fitted-zero coefficients into permanent structural pruning without an update contract. Do not sacrifice equivariant feature flexibility for scalar inference specialization. Keep scalar precontraction optional and general feature/readout parameters intact.
