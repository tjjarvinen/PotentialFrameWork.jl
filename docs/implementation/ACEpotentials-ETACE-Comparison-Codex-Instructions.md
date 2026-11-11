# Fitted Ti–Al–V conversion and CPU benchmark

## Goal

Use the fitted May30 Ti–Al–V model in `data/ref_model` to validate ACEpotentials → ETACE and ACEpotentials → PotentialFrameWork (PFW) conversion. Add an import regression and write a CPU benchmark report for the validated calculators. Stop after the report so the next performance change can be chosen from the results.

Work in the user's local Codex/VS Code checkout. Read `AGENTS.md` and the current implementation status first. Record the starting commit and relevant local changes. Treat the current checkout as authoritative. PR #16 and `benchmark/cpu_ace/DAG-REPORT.md` provide the performance baseline.

## 1. Reload the existing fit

- Locate the training/reconstruction Julia script and fitted `.jls` in `data/ref_model`. Use the actual filenames and checkpoint settings. The expected fit is the reduced degree-10, 500-frame May30 fit, with order 3, cutoff 6 Å, and Ti, Al, and V. Verify these details from the checkpoint.
- Use `load_may30_fit` in a local producer environment pinned to ACEpotentials **0.10.3**. Add only the dependencies needed for the producer and benchmark. Record the actual Julia and resolved package versions, source and checkpoint hashes, and reconstruction commands. Keep the local environment's Project.toml and Manifest.toml out of the committed deliverables.
- Reuse the fitted parameters and state. The original training data is unnecessary for reconstruction and the checks below. Preserve the existing fixture files.
- Reuse the current resolved export workflow where practical. Relevant starting points are `data/ace_ref/export_resolved_model.jl`, `benchmark/cpu_ace/generate_modern.jl`, and `benchmark/gpu_ace/workloads.jl`. The existing recipe-oriented exporter has fixture-specific checks. Factor out the reusable calculator-to-resolved-data logic as needed, preserving its original route.

## 2. Validate polynomial ETACE conversion

Construct the ETACE calculator from the reloaded polynomial model using `ACEpotentials.ETModels.convert2et_full(calc.model, calc.ps, calc.st)`. Retain the one-body and pair components. Verify the effective calculator types, cutoff and species ordering.

Generate fixtures deterministically from a recorded seed or fixed coordinates. Include perturbed bcc, hcp, and ω-like periodic structures, all three species, and an explicitly verified V–V nearest-neighbor case inside the cutoff. Use a common periodic-image regime and check cell face heights and minimum distances. Include an asymmetric displacement and a strain or non-orthogonal cell supported by the compared evaluators. Reuse existing builders and fixture code where suitable. Verify that the generated cases exercise nonzero forces and nontrivial virials.

Evaluate each exact structure using the standard combined energy, forces and virial (E/F/V) interface. Compare polynomial ETACE with polynomial ACEpotentials in Float64. Record maximum absolute energy error per atom, maximum absolute force-component error, and maximum absolute virial-component error, with units, per fixture. Fail on nonfinite values, shape mismatches or agreement outside scale-aware tolerances appropriate for Float64. The previous local five-frame smoke check saw maxima around 2.3e-13 eV/atom, 9.3e-15 eV/Å and 9.8e-14 eV. Use the new comparisons as evidence for the generated fixtures.

## 3. Extend and validate the PFW spline import

Use the reloaded ACEpotentials calculator as the single source for the spline representation. Check the v0.10.3 model-level splinify API, preserving the fitted WB, Wpair and reference energies when constructing the splined calculator. Verify that the splined calculator evaluates with the expected parameter and state layout. Use the actual spline knots and coefficients produced by ACEpotentials, including the radial transforms and envelopes, and record the spline node count.

The model uses **real solid harmonics** and **PolyEnvelope1sR** for the pair basis. Support these in the resolved exporter and PFW importer and evaluator:

- Preserve angular convention, normalization, component order and Cartesian derivatives. Validate solid harmonics against the source basis at several radii and directions so the radius dependence is exercised. Reuse the existing SpheriCart extension and angular dispatch where possible.
- Implement the modern pair envelope and derivative with the same cutoff and near-cutoff behavior as the source. Exercise scalar and SIMD radial evaluation for it. Preserve the separate legacy ACE1 envelope.
- Preserve species, density/product index ordering, coupling, pair weights, reference terms and units. Check that the resolved artifact survives the intended load path. Keep historical resolved files and checkpoint defaults readable. If native calculator saving is supported for the new representation, verify round-trip parity. Inspect shared CPU, functional, precision and device paths for implicit angular or envelope defaults. Any execution path without support for the new representation should fail explicitly.
- Keep exporter and live oracle tooling separate from PFW's core runtime dependencies. Follow existing fixture handling for the resolved artifact and keep its provenance with it.

Add focused regression tests for the new conventions and the normal import/evaluation entry point, including both prefix and DAG E/F/V on the imported model. Use a small deterministic generated fixture set. Add a force finite-difference check and a virial strain check on a representative case, keeping neighbor topology fixed around the perturbations. Run the relevant existing scalar, SIMD, persistence and legacy import tests. Keep the live ACEpotentials/ETACE checks runnable through the dedicated benchmark or integration environment.

Report three comparisons independently:

| Reference | Candidate | What it measures |
|---|---|---|
| Polynomial ACEpotentials | Polynomial ETACE | ETACE conversion parity |
| Polynomial ACEpotentials | Splined ACEpotentials | Spline approximation error |
| Splined ACEpotentials | PFW with the exported tables | PFW import parity |

Choose separate tolerances for numerical parity and spline approximation. If spline approximation is too large for the planned benchmark, inspect radial values and derivatives across the fitted distance range and near the cutoff, then verify an appropriate node count. Preserve the polynomial source and disclose the final node count and errors.

## 4. Benchmark after parity checks pass

- Use CPU Float64 and the same fitted parameters and generated geometry for each workload. Reuse existing benchmark fixtures and harnesses where suitable. Include a small case, an approximately 1,024-atom case and a larger case if runtime permits. Record geometry seed, cell, composition, atom count, cutoff and neighbor counts.
- Measure warmed public combined E/F/V for polynomial ACEpotentials, polynomial ETACE and splined ACEpotentials, and PFW using its imported splines. Record the radial representation for each result. Compare PFW to the splined ACEpotentials path with the same tables. Treat polynomial ETACE timings as a comparison of evaluation paths for the same fit with different radial representations from PFW. Any optional spline ETACE timing needs its own parity check and table provenance.
- Record the exact public-call boundary for each calculator, including graph or neighbor-list construction and reuse, output ownership and allocations. Time calculator reconstruction, conversion, spline generation, plan build and executor preparation outside the warmed call trials.
- Report PFW prepared prefix versus DAG with identical inputs, precision, SIMD and execution settings. Use explicit algorithm settings and record the effective algorithm. Include public PFW timings for both if feasible. Preserve the existing automatic selection.
- Start with one pinned physical core and BLAS threads set to one. If adding a multicore comparison, explicitly verify thread counts and pinning and label serial backends. Use sequential warmed trials, a bounded time budget, preferably ten measured samples per case, and raw samples plus median, spread and allocation data. Report any case that falls short of the sample target.

## 5. Deliver and stop

Write a dedicated Markdown report linked from `benchmark/cpu_ace/README.md`, and update the implementation status with verified coverage and any gaps. Store machine-readable timing and numerical records under a new result directory alongside exact reproduction commands. Preserve earlier benchmark reports.

Include model provenance and actual basis/product counts, the three numerical comparisons, public E/F/V medians and timing boundaries, prepared prefix versus DAG, and the relevant regression test outcomes. Explicitly identify missing or failed gates and omit performance ratios for calculators that failed parity. Keep this milestone scoped to ACEpotentials, ETACE and PFW. Stop at the validated baseline and report.
