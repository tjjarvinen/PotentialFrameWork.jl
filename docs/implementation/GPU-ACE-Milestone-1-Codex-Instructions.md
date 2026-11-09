# GPU ACE milestone 1 — Codex implementation instructions

## Objective and authority

Implement scalar ACE energy, forces, and virial on AMD GPUs in PotentialFrameWork.jl. Add model precision/device adaptation, a precision-comparison utility, a reusable benchmark suite, and a CPU SIMD audit report.

These instructions record the user's agreed scope. They authorize GPU ACE and supersede older phase-specific prohibitions on GPU ACE, required GPU dependencies, and GPU atomic force reduction where necessary for this milestone. Preserve all unrelated repository instructions, package identity, supported CPU behavior, units, parameter authority, output ownership, and existing CPU AD contracts. Keep AtomsCalculators and the local neighbor protocol as interfaces. Execution and allocation belong in implementation modules.

Repository source is https://github.com/tjjarvinen/PotentialFrameWork.jl. Planning was grounded in commit `5eef1c1` (General equivariant readout for ACE). Inspect the actual local branch and uncommitted files before editing. Build on subsequent changes if present. Preserve user-provided data and pinned reference environments. Read the current AGENTS.md, relevant implementation specifications, and STATUS.md before implementation.

Complete routine reversible work autonomously. Resolve implementation details through source inspection and measurement. If hardware, model compatibility, or missing inputs block a gate, complete independent work and report the precise blocker. Record actual test results and limitations. Do not claim AMD validation from a CPU execution of a portable kernel.

## Milestone deliverables

- GPU scalar ACE evaluation for every supported nonempty subset of Energy, Forces, and Virial, with explicitly defined behavior for an empty request.
- Float64 and Float32 numerical model data, parameters, computation, and buffers, with explicit conversion and backend selection.
- Optimized imported spline radials and custom radial callables.
- A general staged GPU path, a measured compact fused candidate, and a one-thread-per-site benchmark baseline.
- Minimal site intermediate storage, optional site batching, and atomic force accumulation.
- AMD hardware correctness and performance results.
- CPU/GPU precision-comparison utility accepting user structures.
- Quick/full BenchmarkTools.jl presets using the existing TiAl model and public Ti and Ti–Al–V models, including independent ACEpotentials comparisons.
- CPU SIMD audit and report using the same benchmark workloads where feasible.

Equivariant GPU outputs, GPU training, GPU AD, force-loss derivatives, full MACE layers, MPI, long-range electrostatics, multiple periodic images on GPU, and NVIDIA/Apple hardware validation are later milestones. Maintain architectural boundaries suitable for future pullbacks and feature readouts. Preserve already supported CPU equivariant and AD functionality.

## Source orientation

Inspect these components and their current callers before changing representations.

- `src/submodules/ACE/cpu.jl`, `serial.jl`, `projection.jl`, `resolved.jl`, `radials.jl`, `calculators.jl`, and persistence/interface files.
- `src/submodules/ACE/features.jl` and `feature_coupling.jl` for future reuse boundaries and CPU regression coverage.
- `src/submodules/NeighbourList/gpu.jl`, `execution.jl`, and the local `NeighbourInterface`.
- `src/submodules/PairPotentials/gpu.jl` for existing GPU preparation, backend, result, and buffer conventions.
- `ext/PotentialFrameWorkSpheriCartExt.jl` and `src/submodules/ACE/angular_blocks.jl`.
- `data/ace_ref`, particularly `export_resolved_model.jl`, saved fixtures, and independent validation scripts.
- `test/gpu_pairs.jl`, existing ACE/AD/interface tests, `benchmark/ace`, and GPU pair/skin benchmarks.

At the planning baseline, the existing resolved exporter is deliberately pinned to one ACEpotentials 0.10.2 model, with assertions for 270 parameters. It is not a general importer. The CPU executor stores Float32 scratch with Float64 imported coefficients in its current mixed-precision mode. Both facts need explicit attention in this work.

## 1. Model preparation, precision, and storage adaptation

Expose an array/storage constructor keyword and a separate precision keyword through the existing relevant constructors/import paths. Illustrative usage is `array_type=Array, precision=Float64` for CPU and `array_type=ROCArray, precision=Float64` for AMD. Exact names should follow repository conventions. CUDA and Apple array types describe future interface compatibility, not validated support in this milestone.

Provide an explicit operation that creates an independent adapted model from an existing model, selecting array type and precision. Use the same adaptation machinery for construction and conversion. Reuse Julia's existing adaptation mechanisms where suitable. Avoid a parallel family of nearly identical public constructors.

- Convert floating-point coefficients, spline knots/tables, transforms, envelopes, numerical coupling weights, references, explicit parameters, and custom configuration where appropriate.
- Preserve species identities, index arrays, shapes, masks, product structure, and normalization/order metadata. Distinguish integer numerical metadata from floating-point coefficients.
- Move kernel-accessed numerical tables and indices to the selected device. Descriptive and construction metadata may remain on the host.
- Prepare GPU-safe packed representations of sparse/ragged tables where host representations cannot execute on device. Retain a single authoritative mathematical model and parameter ordering.
- Infer backend from actual numerical array storage, including the selected device. Verify compatible model, parameters, positions, output, and scratch placement. Make adaptation/transfers explicit at preparation or input-update boundaries.
- Repeated execution must not silently rebuild or transfer the model. Parameter changes must reach kernels without stale precontracted weights.
- Keep existing explicit parameters/state semantics. If parameters live outside the model, adapt the complete model/parameter bundle coherently.
- Preserve ownership across source/adapted models and independently prepared executors. Numerical results that are documented as owning must survive later calls and buffer growth.
- Conversion is an explicit preparation operation. Preserve existing CPU functional AD entry points and avoid inserting conversions that detach their inputs.
- Document that Float32-to-Float64 conversion cannot recover discarded coefficient precision.

Preserve core Julia compatibility where practical. Put device-specific integrations in optional dependencies/extensions or the existing backend integration structure. Keep loading and using CPU functionality possible without AMDGPU installed.

## 2. GPU radial and angular evaluation

### Imported splines

Spline evaluation is the primary production optimization target. Preserve the actual imported piecewise cubic function, transform, envelope, cutoff, and analytic derivative. Do not substitute the generating analytic polynomial for an imported spline. Convert data to the selected precision consistently and validate interval boundaries and envelope derivatives.

### Custom radials

Support custom radial callables through the existing radial protocol where practical. They return/write a fixed declared channel count and first derivatives with respect to distance, using explicit parameters/configuration and known scratch requirements.

Branching is allowed. Treat numerical smoothness, expensive evaluation, divergent branches, unsupported allocations, and unsupported host operations primarily as documented user responsibilities. The evaluator must still enforce structural requirements needed for correct execution. Do not build a general callable analyzer, blacklist branching, or require finite-difference sampling at every preparation. Optional development checks may reuse existing radial diagnostics.

Allow CPU-only custom callables to continue working on CPU. For GPU, provide clear errors for compilation/adaptation failures and document any user extension hook for custom configuration types. Automatic differentiation of arbitrary custom radial functions is deferred. Radials that supply explicit first derivatives are sufficient here.

Protect inactive/padded candidates before operations such as reciprocal distance or angular normalization. Apply actual physical cutoffs even when cached connectivity contains skin candidates. Use an explicit coincident-atom policy compatible with existing CPU behavior.

### Angular evaluation

Use a verified device-compatible angular implementation. Check the actual installed SpheriCart APIs before attempting reuse. Preserve real-harmonic normalization, component ordering, derivatives, and the lmax=0 special case. Reuse existing recurrence logic where correct and useful. Prove angular values and Cartesian derivatives against the independent CPU reference across supported lmax and precision.

## 3. Geometry and neighbor integration

Reuse existing GPU sorted cells, permutation/species mapping, device buffers, directed connectivity, and optional skin. Preserve orthorhombic/single-image restrictions and existing geometry validation. Topology and periodic geometry remain fixed during an evaluation.

ACE needs complete directed neighborhoods. Preserve central selection, species-dependent channels, supported descriptor masks/pruning and species/order filtering, physical cutoffs, fixed references, and the imported pair component. Confirm every current supported scalar-model feature has a correct GPU implementation or a precise documented limitation. Avoid silently dropping terms or copying the calculation to CPU.

Validate sorting and restoration to original atom order, species/metadata changes, cell changes, neighbor rebuilds, buffer growth, empty systems, isolated atoms, and shrinking workloads. Reuse pair/neighbor infrastructure without regressing its behavior.

## 4. Thread mapping and partial fusion

Implement the general GPU computation in three main mathematical stages. Auxiliary initialization and final reductions may require additional kernels.

1. Evaluate radial/angular quantities and accumulate site densities. Start with a workgroup per site, using bounded neighbor/channel tiles and cooperative reductions.
2. Evaluate site products and scalar readout, accumulating energy and density adjoints in the same stage. Tile products/channels as necessary. Preserve the scalar coefficient-precontraction optimization where mathematically applicable.
3. Revisit edges, recompute radial/angular values and first derivatives, contract density adjoints into edge gradients, and accumulate forces/virial.

The stages may choose different thread mappings. A large basis must not require all channels to fit in shared memory. Define ownership and synchronization explicitly. If several workgroups contribute to one site, complete the density before contraction through valid global staging/kernel boundaries. Avoid races between product adjoint updates and between atom force updates.

Add a compact candidate using one cooperative subgroup per site and several sites per workgroup where appropriate. Fuse the site calculation when local state permits. Use backend-portable workgroup primitives or verified backend abstractions rather than assuming NVIDIA warp width or CUDA-only intrinsics.

Keep one-thread-per-site evaluation as a correctness/performance baseline. Reuse arithmetic and packed tables across strategies. Avoid maintaining unrelated mathematical implementations for each mapping.

Provide a manual strategy override and simple benchmark-derived automatic thresholds selected during executor preparation. Candidate names such as `:auto`, `:staged`, and `:compact` are illustrative. Use model dimensions, precision, backend, and only measured workload factors. Sophisticated autotuning and per-step strategy switching are unnecessary.

Benchmark the compact candidate. Promote it to automatic selection only where measured benefits justify it. If it is consistently slower, retain a documented experimental override only if useful, or remove it and record the result. A validated staged path plus a documented negative compact experiment satisfies this aspect of the milestone.

## 5. Intermediate storage and reduction

Default to reusable site densities and their adjoints, with bounded product scratch. Preserve extra density components needed by existing species/order filtering. Recompute edge basis quantities in the force pass. Avoid global arrays of every radial/angular value and derivative per edge in the default path.

Keep density and adjoint buffers separate initially. Report memory using actual counts of density components, products, and policy-specific buffers. Distinguish these from radial/output channels and final invariant basis counts.

Support optional batches of central sites to bound scratch memory. Batches complete the relevant three stages and accumulate into common final outputs. Batch boundaries must not lose contributions to atoms outside the central batch. Connectivity, final outputs, and model data remain additional memory costs beyond site scratch.

Use atomic force accumulation. Reduce each central site's contribution locally before atomics where practical. Handle neighbor scatter correctly. The gradients associated with directed edges i→j and j→i derive from different site energies and generally differ. Do not infer one from the other by a sign flip.

Use per-site/workgroup partial sums and final reductions for energy and virial where suitable. Preserve existing force signs, virial conventions, pair counting, and units. No per-worker full-system force replicas. Edge-gradient storage/gathering is deferred pending profiling.

Energy-only evaluation should avoid adjoint and edge-force work. Virial-only evaluation still requires the relevant displacement derivatives. Reset every requested output correctly on repeated calls. Measure reduction-order differences with appropriate tolerances rather than requiring bitwise equality across mappings.

## 6. Public benchmark models and provenance

Use the existing small TiAl reference plus these public fitted models.

| Model | Reported invariant basis count | Published settings |
| --- | ---: | --- |
| Pure Ti | 1,809 | Correlation order 3, total degree 20, cutoff 6 Å |
| Ti–Al–V | 29,703 | Correlation order 3, total degree 15, cutoff 6 Å |

Source archive is https://zenodo.org/records/15114121, DOI 10.5281/zenodo.15114121. Paper is Allen and Bartók, “Multi-phase dataset for bulk Ti and the Ti-6Al-4V alloy”, Scientific Data 12, 1925 (2025), https://doi.org/10.1038/s41597-025-05302-3.

Relevant files include

- `Ti_combined_data_oct24_2023_ace_order_3_degree_20_rcut_6.0_cw_moreV.json`
- `ti64_combined_data_Apr08_rcut_6.0_order_3_degree_15_custom_reg.json`
- Corresponding extended-XYZ data and YACE/pair-table exports.

The May30 alloy model is an alternative fit with the same published basis settings. Choose one alloy fit, record it, and pair it with the appropriate structures. Do not download the unrelated GAP model files.

The archive API metadata was checked during planning and declared CC BY 4.0. Recheck the record and downloaded notices during acquisition. Retain authorship, source DOI, British Crown/AWE copyright notices, license/link, file checksums, and a description of conversion changes. Keep third-party artifact licensing separate from package code licensing. Provide deterministic download/preparation commands and cached, hash-verified data. Prefer explicit benchmark-data acquisition over downloads during normal package loading or core tests.

### Compatibility gate

Do not assume these JSON files use the current TiAl schema. Inspect source model format, original package versions, basis definitions, radials, pair terms, references, and normalization first. Establish an isolated reproducible producer environment and independent ACEpotentials oracle.

Implement a narrowly scoped faithful exporter/adapter if needed. Preserve original basis and fitted parameters. Extend the resolved representation only where required and keep existing checkpoint/import compatibility. Avoid refitting models, substituting random coefficients, or comparing different physical functions under the same benchmark label.

If the models require a new radial family, explicitly implement or faithfully tabulate it with measured approximation error. Separate conversion approximation from Float32 rounding and GPU implementation errors. If import is blocked, report exactly which feature/version prevents it and continue independent GPU work with TiAl. Do not declare the public-model benchmark gate complete until actual model predictions are verified.

## 7. Precision-comparison utility

Provide a reusable function accepting one model/parameter bundle and either one structure or several user-supplied structures. Array type/backend and requested properties are configurable. Return a structured report with a readable display. Public names should follow existing conventions.

Create Float64 and Float32 versions from the same source model. Prefer original Float64 coefficients. Preserve and report the provenance precision if only a Float32 source exists. Construct both geometries from the same source coordinates/cell and report rounding to selected precision.

Report per structure and aggregate over structures

- Absolute total-energy difference and per-atom energy difference.
- Force RMS component error, maximum absolute component error, and maximum per-atom vector difference.
- Virial RMS component error and maximum absolute component error. Optional stress comparison requires a well-defined nonzero volume and explicit units/convention.
- Nonfinite values and evaluation failures, associated with structure and precision.
- Model/structure identity, evaluated atom counts, precision/backend, requested quantities, tolerance settings, and worst-case structure/component identifiers.

Use absolute errors as primary metrics. Define any relative normalization explicitly and handle near-zero references without misleading ratios. Handle zero atoms, absent properties, isolated atoms, and empty structure collections deliberately. Document aggregate weighting, distinguishing pooled component RMS from averages of per-structure RMS.

With no tolerances, provide descriptive numerical status and differences. With explicit tolerances, add property/structure pass/fail results. Nonfinite outputs and evaluation errors must be visible failures and must never improve aggregate results through silent exclusion. Keep failed cases identifiable if processing continues.

The utility measures precision sensitivity on supplied examples. It does not establish universal Float32 suitability. Benchmark presets run numerical checks before timing and carry descriptive precision errors into their output. A default Float32 suitability threshold is unnecessary. Mandatory implementation-correctness checks remain separate from user-selected application tolerances.

## 8. BenchmarkTools suite and user workflow

Add BenchmarkTools.jl in a dedicated benchmark environment, reusing existing benchmark scripts/data where practical. Keep heavy oracle/producer environments separate if dependency versions conflict. Users should be able to run the suite on their computer without editing source code.

Provide documented commands for data preparation and quick/full runs, backend/model selection, precision, strategy override, sample/time budget, and output destination. Avoid a combinatorial default matrix. Proposed case sizes below are starting defaults to verify for runtime and memory, not previously measured thresholds.

### Quick preset

- Cover all three models with a bounded set of deterministic structures, including at least one 10k+ atom workload when resources permit.
- Measure prepared E/F/V, a representative complete calculator call, and the available independent CPU oracle.
- Include CPU and AMD when available. CPU-only invocation must work and clearly identify omitted GPU cases.
- Use brief explicit timing budgets and record actual elapsed runtime. Avoid rerunning identical expensive oracle cases redundantly.

### Full preset

- Extend atom counts, neighbor density, one/multiple species, Float32/Float64, CPU SIMD widths/thread counts, and staged/compact/baseline GPU comparisons.
- Include prepared evaluation, complete application calls, and a deterministic sequence of displaced structures with cached-neighbor reuse and forced rebuilds.
- Record workload and stage-level diagnostics supporting dispatch thresholds and the SIMD audit.

Build large orthorhombic structures from source data or physically appropriate bulk cells with deterministic perturbations. Include small independent-oracle fixtures, all species, and varying neighborhoods. Replicate cells enough to satisfy the GPU single-image constraint. AtomsBuilder/AtomsSystems may be used in the benchmark environment. Avoid including structure-building time in prepared evaluation.

### Timing and comparison contracts

- Compare the same physical models on identical structures. Verify E/F/V before reporting performance ratios.
- Warm each configuration before timing. Record setup/compilation separately if measured.
- Synchronize GPU work before timing and ensure the timed operation completes before the sample ends. Do not report launch latency as completed evaluation time.
- Use `evals=1` for stateful/rebuild cases and control resets in sample setup. A sequence benchmark must replay the same sequence for every sample, not drift continuously across samples.
- Document the exact timed boundary. Prepared execution keeps model/inputs/buffers resident. Complete calls include transfers/preparation actually performed by that public call. Report transfers explicitly.
- Use a common Float64 reference for fair correctness comparisons. Label mixed-precision paths and unsupported oracle precision combinations accurately.
- Run CPU thread variants in separate Julia processes, record BLAS/Julia threads, avoid oversubscription, and run competing configurations sequentially.
- Where oracle environments conflict, use separate processes with the same serialized fixtures and model provenance.
- Report median, minimum, spread/sample count, time per evaluation, and atom throughput. Retain raw samples or BenchmarkTools trials where practical.
- Record CPU, GPU, driver/runtime, OS, Julia/package versions, repository commit/dirty state, model/structure hashes, basis and density sizes, neighbor counts, skin/rebuild behavior, precision, strategy/tile settings, threading, and timing configuration.
- Distinguish host allocation measurements from device buffer sizes or backend-specific device-allocation measurements. Keep unavailable metrics marked unavailable.
- Output a readable summary and a versioned machine-readable report, using existing TOML conventions or a similarly simple format.

## 9. CPU SIMD audit

This deliverable is an investigation and report. Plan subsequent CPU optimization from its findings. Limit source changes to necessary benchmark/inspection hooks or correctness fixes.

Inspect explicit SIMD and compiler-generated vectorization separately for geometry, spline transforms/envelopes/lookup, spline Horner, angular values/derivatives, density projection, products/reverse, edge-gradient contraction, and force/virial reduction. Include custom-radial lane-wise behavior and the separate CPU equivariant feature path as clearly identified observations where relevant.

At the baseline, imported cubic Horner uses explicit SIMD, transforms/interval selection are scalar per lane, custom radials are evaluated lane-wise, and optional angular blocks/fusion already exist. Do not equate scalar-looking source with absence of compiler vectorization. Existing batched SpheriCart has beaten the experimental explicit blocks in recorded cases.

Inspect optimized LLVM and native code for representative warmed concrete signatures. Identify packed arithmetic, scalar arithmetic, gathers, reductions, allocation, and relevant layout constraints. Record hardware target, Julia version, model, precision, and width. Static evidence alone does not establish runtime SIMD percentage.

Measure scalar/SIMD modes, existing angular/fusion alternatives, both precisions, and representative models/neighborhoods. Report full E/F/V timings alongside kernel diagnostics, including species packing, tails, and phase occupancy where they affect conclusions. Clearly label allocations and preparation cost. Suggest measured priorities such as density/force contraction or SIMD across central sites. Preserve all original measurement artifacts and label new ones by environment and commit.

## 10. Implementation sequence

1. Inspect current state and create a short repository-local implementation record. Resolve model formats and provenance while keeping existing fixtures intact.
2. Implement coherent precision/storage adaptation and GPU-safe model tables. Validate CPU round trips and independent ownership first.
3. Implement GPU radial/angular primitives and establish AMD device tests.
4. Implement the staged scalar E/F/V path with correct directed geometry, filters, pair/reference terms, atomics, and property requests.
5. Add reusable storage, optional site batching, repeated-call/growth handling, and resident input/output workflow.
6. Implement the precision-comparison utility and public-model oracle fixtures/export support.
7. Build quick/full benchmarks and evaluate baseline, staged, and compact candidates. Establish simple selection thresholds only from measurements.
8. Complete the CPU SIMD audit and documentation. Run targeted regressions and finish the status/validation record.

Stages may overlap when independent. Keep numerical correctness gates ahead of speed claims.

## 11. Validation and completion gates

### Numerical behavior

- Float64 AMD E/F/V agrees with independent CPU evaluation for TiAl, Ti, and Ti–Al–V across representative structures. Validate conversion separately before using a converted model as the GPU oracle.
- Check force components against energy-position finite differences and virial against homogeneous-strain energy derivatives, maintaining the project's conventions and safe fixed-topology perturbations.
- Check translation invariance, energy rotation invariance, force/virial covariance, atom permutation, species mappings, and supported selection/filtering behavior.
- Cover cutoff boundaries, spline knots, custom analytic radials and their derivatives, lmax=0 and higher lmax, empty/isolated cases, and physical zero-distance handling.
- Use property-specific absolute/relative tolerances motivated by scale and precision. Record chosen tolerances and worst errors. Account for fixed reference offsets and cancellation. Do not loosen a tolerance merely to hide a failure.
- Float32 AMD implements the correctly rounded model consistently within justified implementation tolerances. Report its application-level differences from original Float64 separately.

### State and storage behavior

- Repeated evaluation resets and accumulates correctly, handles neighbor reuse/rebuilds and changed parameters, and restores original order.
- Sorting, skin, species changes, cell changes, growing/shrinking capacity, all supported property subsets, and multiple site batches are covered.
- Forced staged/compact/baseline results agree within precision-aware reduction tolerances. No silent CPU fallback in device tests.
- Verify model adaptation independence, source immutability, all floating tables' types, and output lifetimes. Preserve supported checkpoint loading and CPU AD regression tests.
- Report resident model, connectivity, outputs, and scratch storage separately. Warm prepared execution reuses allocations within capacity. Any unavoidable allocation or synchronization is documented.

### Usability and evidence

- AMD hardware tests run and record hardware/software versions. CPU-only load/tests work without AMDGPU. Existing pair GPU and CPU functionality passes relevant regressions.
- Both presets have reproducible commands. Timing boundaries and exclusions are explicit. Public-model acquisition/export is documented with hashes and attribution.
- Precision reports accept a single structure and several structures, handle failure/nonfinite cases, and provide descriptive and tolerance-based modes.
- The CPU SIMD audit report includes source and generated-code evidence, measurements, limits, and follow-up priorities.
- No universal speedup or mandatory performance ratio is required. Automatic compact selection requires demonstrated benefit. Report regressions and negative experiments honestly.
- Update relevant docs and STATUS.md with completed gates, exact reproduction commands, remaining limitations, and measured results. Include third-party notices and report schemas.

If required AMD hardware is unavailable in the implementation environment, supply runnable AMD tests and complete independent work, but mark hardware gates unverified. Likewise, a model-import blocker leaves the corresponding benchmark gate incomplete. Do not substitute a successful synthetic or CPU-only case for the required evidence.

## Reference links

- PotentialFrameWork.jl — https://github.com/tjjarvinen/PotentialFrameWork.jl
- Model archive — https://zenodo.org/records/15114121
- Model paper — https://doi.org/10.1038/s41597-025-05302-3
- CC BY 4.0 terms — https://creativecommons.org/licenses/by/4.0/
- BenchmarkTools manual — https://juliaci.github.io/BenchmarkTools.jl/stable/manual/
- KernelAbstractions API — https://juliagpu.github.io/KernelAbstractions.jl/stable/api/

Check actual installed dependency versions and APIs during implementation. This document specifies the agreed implementation work, not completed code or benchmark results.
