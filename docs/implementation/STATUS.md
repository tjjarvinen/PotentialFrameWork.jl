# Implementation status

## CPU ACE interface phase (2026-10-24)

Started on `ACE_tuning` at `71492894e2e3881c75bc032dacab2879a4df0d57`, identical
to local `main`; preserved the user's modified AGENTS.md and untracked interface
specification/plan/validation documents. Package name/UUID/authors/license and
Julia 1.10.10 compatibility are unchanged. Source remains in ACE, the root stays
thin, and no GPU dependencies or implementation were added.

The public API and mathematical conventions are documented in
`docs/ace_interface.md`; runnable examples are `examples/ace_interface.jl`,
`examples/ace_radial_source.jl` and the separate pinned-environment
`examples/acepotentials_radial_adapter.jl`.

| Milestone | Implementation and evidence |
| --- | --- |
| I0 | Inventory of ResolvedModel dispatch in numerical/calculator/CPU/ChainRules/Lux/fitting/checkpoint code; HEAD equals main; baseline 44,501 assertions on Julia 1.13.2 / 4 threads (`/tmp/pfw-interface-baseline.log`). Reference hashes below were preserved. |
| I1 | Owning `descriptors`/`basis_summary`, namespaced original IDs and actual mixed-order coupling terms; immutable `BasisSelection`; reversible effective masks; structural pruning with coefficient/dependency maps and nonlinear scale transfer; `restrict_parameters` supports reverse scatter into parent coordinates. `restrict_basis(path; ...)` compacts serialized radial/coupling/parameter tables before constructing ResolvedModel, without constructing a full evaluator/workspace. Mask/prune/update/empty-basis, mapped prior, mixed-order metadata and repeated-provenance tests pass. |
| I2 | Original-order `centers` at numerical, calculator, prepared CPU and fitting entry points; Lux inputs carry centers and ChainRules owns saved selections. All/range/reordered/empty/invalid behavior, site partition sums, forces on nonselected neighbors, force-write storage and saved selections are checked. Full environments are still prepared. |
| I3 | One selected neighbor-species group and inclusive factor-order range via `project_species`; U+tV truncated convolution and division-free reverse, shared readouts, explicit reference/pair order accounting and reusable per-slot polynomial scratch. Independent tiny-factor enumeration, zeros/repeated factors, sum-of-orders, full-minus-excluded, rotation, strain, combined mask/center/design/AD checks pass. |
| I4 | Concrete radial channel families, explicit `RadialMetadata`, in-place `CallableRadial`, analytic/value-only and allocating adapters, separate custom ps/config, declared capability diagnostics and local ForwardDiff second/mixed probes. Imported SIMD splines retain their specialized path. Smooth trainable custom examples pass scalar/Vec Float32/64/mixed-precision comparisons and force-loss position/parameter products. Pinned ACEpotentials 0.10.2 `Models.evaluate_ed` adapter passes 35 checks (`/tmp/pfw-interface-pinned-adapter.log`). |
| I5 | Native checkpoint v2 stores policies/provenance/custom definitions and authoritative ps/st once; old v1 native/resolved loading remains tested. Registered factories, explicit trusted-source bundles and opt-in raw-distance Hermite snapshots round-trip. Source refusal precedes execution, namespaces are fresh, world transitions stay at public execution/loading boundaries, missing dependencies/hash/entrypoint/interface/reference failures are explicit. Fresh-process relocated source loading, repeated custom-pruning checkpoints and a code-disabled sentinel test pass; invalid saves preserve existing files. |
| I6 | UnicodePlots 3.9 weakdep/extension using its direct API (no Plots.jl): order/index counts, coefficient distribution, radial curves, coupling sparsity and grouped counts. Six plot types including empty bases pass 17 checks; actual ACEfit 0.3.1 transformed Unitful-label prediction recovery and non-diagonal rectangular mapped-prior solves add four checks (`/tmp/pfw-interface-optional.log`). Core-only transformed CPU/checkpoint round trip passes with plotting, ACEpotentials, SpheriCart and AD weak dependencies absent (`/tmp/pfw-interface-core-only.log`). |
| I7 | Full supported-version integration, small allocation/size/overhead measurements, examples and final reproduction evidence are recorded below. |

Numerical, persistence and AD boundaries retain full-system original-order owning
outputs, all directed periodic images, imported pair-site accounting, fixed
cell/topology differentiation, saved-input ownership and the existing optional
virial-position products. Production uses one sorted force array and exclusive
scratch slots; masks/readout caches and center flags are refreshed in reusable
executor buffers. No model-derived precontraction is authoritative.

Material limits: callable CPU execution currently uses scalar lanes; direct
SIMD-compatible callable inputs are checked separately. Scalar AD/allocating
adapters and pruned custom-channel gathers can allocate and recompute parent
channels. Noncopyable mutable callable captures require explicit config/ps;
immutable code/module references are retained while numerical state is copied. Angular backends still compute all harmonics through the original
lmax. Custom parameter gradients use local chunked ForwardDiff in analytic
site reverse products; no whole-system Hessian is materialized. Capability probes
are sampled checks, not proofs for arbitrary callables/backends/AD orders.
Registered channels must be available at load; source loading requires explicit
`allow_code=true`, is not a security sandbox and never installs dependencies.
Tabulation is an explicit Float64 fixed-parameter approximation with off-grid
value/derivative errors; it is C1, generally not C2 at knots, and does not preserve
original custom-parameter gradients. Unknown custom index weights require an
explicit prior or the explicit neutral prior option. General basis generation,
independent species counters/coefficients, changing regions, site-energy APIs,
composition/QM-MM physics, cell AD and GPU/MACE remain outside this phase.

Final complete suites pass **44,757 assertions per process** on **Julia 1.13.2 /
4 threads** and **Julia 1.10.12 / 4 threads**. Logs:
`/tmp/pfw-interface-final-113-t4.log`, `/tmp/pfw-interface-final-110.log`.
The standalone final interface suite passes **255 checks on Julia 1.13.2 /
1 thread** (`/tmp/pfw-interface-final-113-t1-focused.log`). This includes ordinary
and scalar-AD radial force losses, nested custom-pruning checkpoint loads and
rebound custom channel compatibility. No required interface milestone remains
unimplemented; the explicitly deferred/optional capability limits above remain.

The pinned reference revalidation passes **1,276 checks**: 1,136 radial
value/derivative comparisons, 132 checks on all 33 original-evaluator E/F/W
predictions and eight original-coordinate design/ACEfit checks
(`/tmp/pfw-interface-reference.log`). Maximum absolute discrepancies are
**1.4551915228366852e-11 eV**, **6.146194664324867e-13 eV/Å** and
**8.43769498715119e-13 eV**. The source model was loaded, never overwritten or
refitted. The optional scalar interoperability adapter separately passes 35
checks. The optional plotting/fitting suite passes **21 checks per process** on Julia
1.13.2 and 1.10.12 (`/tmp/pfw-interface-optional.log`,
`/tmp/pfw-interface-optional-110.log`),
including actual QR solves with non-diagonal rectangular `R_old*T` operators;
ACEfit 0.3.1's square-prior target-length assumption is handled by augmenting the
objective explicitly in our extension, then calling its matrix QR solve.

Hardware: AMD Ryzen 9 9950X, 16 cores / 32 logical CPUs. Versions: AtomsBase 0.5.3,
AtomsCalculators 0.2.3, SIMD 3.8.1, Unitful 1.29.0, SpheriCart 0.2.4,
ForwardDiff 1.4.6, ChainRulesCore 1.26.4, Zygote 0.7.14, LuxCore 1.5.3;
UnicodePlots 3.9.0 and ACEfit 0.3.1 in the optional application environment.
Core Julia 1.13 tests use StaticArrays 1.9.23 / NeighbourLists 0.6.2;
Julia 1.10 uses StaticArrays 1.9.24 / NeighbourLists 0.5.10. The application and
pinned reference environments use StaticArrays 1.9.24. Julia 1.10 source/tests
were copied to `/tmp/pfw-interface-110`, retaining its separately resolved
compatible Manifest. Root Pkg.test emits a stale ignored-Manifest project-hash
warning after the weakdep addition; the test sandbox resolves successfully. The
pinned reference Project/Manifest were not resolved or changed.

Reproduction (application environment `/tmp/pfw-interface-env` contains
PotentialFrameWork as a path dependency plus AtomsBase, Unitful, SIMD,
SpheriCart, ForwardDiff, ChainRulesCore, Zygote, LuxCore, ACEfit and UnicodePlots):

```sh
julia --project -e 'using Pkg; Pkg.test(; julia_args=["--threads=4"])'
julia +1.10 --project=/tmp/pfw-interface-110 -e 'using Pkg; Pkg.test(; julia_args=["--threads=4"])'
julia --project=/tmp/pfw-interface-env --threads=1 test/interface_runtests.jl
julia --project=/tmp/pfw-interface-env --threads=4 -e 'using PotentialFrameWork, SpheriCart, Test; const ACE_MODEL_PATH=joinpath(pwd(),"data/ace_ref/resolved_model_v1.toml"); include("test/ace_interface_optional.jl")'
julia +1.10 --project=/tmp/pfw-interface-optional-110 --threads=4 -e 'using PotentialFrameWork, SpheriCart, Test; const ACE_MODEL_PATH=joinpath(pwd(),"data/ace_ref/resolved_model_v1.toml"); include("test/ace_interface_optional.jl")'
julia --project=data/ace_ref --threads=4 data/ace_ref/validate_resolved_model.jl
julia --project=data/ace_ref --threads=4 examples/acepotentials_radial_adapter.jl
julia --project=/tmp/pfw-interface-env --threads=4 examples/ace_interface.jl
julia --project=/tmp/pfw-interface-env --threads=4 benchmark/ace/interface.jl
```

The small diagnostic artifact is `benchmark/ace/interface_measurements.toml`;
its reproduction script and workload are documented in `benchmark/ace/README.md`.
On three periodic atoms / 56 directed active edges, serial prepared E/F/W with
four exclusive scratch slots and one BLAS thread records **96 bytes allocated
per prepared execution** for all five policies (not a zero-allocation claim).
Seven samples of 100 repetitions report median times of **36.34 μs** full,
**35.51 μs** masked, **17.63 μs** pruned, **32.74 μs** order 0 and **36.95 μs**
orders 1:2; min/max variation is retained in the artifact. The pruned basis has
20/20/12/12 invariant/product/density/radial dimensions versus 123/167/27/30;
both still compute nine harmonics. Model sizes are 560,904 / 316,240 bytes
(full/pruned), with 108,360 / 73,160 bytes of per-slot scratch summed over four
slots. The sole production force matrix is **72 bytes** in every case, independent
of the four slots; caches, parameter trees and total Julia executor sizes are
reported separately. Owning numerical outputs allocate 336 bytes;
full boundary preparation allocates 56,880 bytes (56,432 pruned). These numbers
are a small geometry-specific diagnostic, not a general speedup or end-to-end MD
performance claim. Environment construction and owning outputs are excluded from
the prepared kernel times and reported separately as allocation costs. Initial
per-call mask/readout temporaries were replaced with owned reusable buffers before
these final measurements. `git diff --check` passes.

Reference SHA256 values (unchanged):
- model: `a61a9e57c91cc25663359dd876c47fe9a4c3662ee8d5b4a27a92350c90435f0e`
- XYZ: `249c00cc8ce9a2ca6e9d45f74b99c9690a4ade5d00df976deb4fd1365737c62c`
- reference Manifest: `0aeb69de6df1ce209afcc44eec338a3dd632db02801a7902324945fbbb8b7ed2`
- reference Project: `be49f9b97d798ad7c66512e8abc7af42d729fc41a47cb7a3b3f4d72f1d06d68c`

## Unitful ACE training labels (2026-10-23)

`ACE.training_rows` now converts Unitful energy and configurational virial labels
to eV and force labels to eV/Å before subtracting fixed baselines. Plain numbers
retain the same eV/Å conventions. Compatible scaled units (including J, N, meV
and eV/nm) and mixed plain/unitful properties are accepted; incompatible
dimensions raise Unitful.DimensionError. Residual weights retain their existing
meaning in the converted units. Force labels also accept vectors of per-atom
three-component vectors, preserving original atom/component ordering for
ExtXYZ/AtomsCalculators outputs. Missing labels and `nothing` still omit rows.
The AtomsBase geometry boundary already converts lengths to Å; explicit numerical
kernels/design matrices continue to use plain arrays. No dependencies changed.
README and the training example document and demonstrate Unitful labels.

Validation adds 38 checks for conversions, component ordering, mixed inputs,
dimension/length errors, omitted rows, Float32/empty force labels, input ownership,
and design-matrix prediction recovery with nm AtomsBase geometry and Unitful
calculator labels for both jointly fitted and fixed pair terms. The pre-edit
Julia 1.13.2 / 4-thread suite passed 44,463 assertions
(`/tmp/pfw-unitful-baseline.log`); the post-edit suite passed 44,501
on both Julia 1.13.2 and 1.10.12 with four threads
(`/tmp/pfw-unitful-final-113.log`, `/tmp/pfw-unitful-final-110.log`). Hardware:
AMD Ryzen 9 9950X, 16 cores/32 logical
CPUs. Unitful 1.29.0, AtomsBase 0.5.3, AtomsCalculators 0.2.3, SpheriCart 0.2.4,
SIMD 3.8.1, ForwardDiff 1.4.6, Zygote 0.7.14. Reproduce with
`julia --project -e 'using Pkg; Pkg.test(; julia_args=["--threads=4"])'`.
Julia 1.10 used an isolated copy at `/tmp/pfw-unitful-110` with its own resolved
Manifest (StaticArrays 1.9.24, NeighbourLists 0.5.10; Julia 1.13 uses 1.9.23 and
0.6.2). Reproduce with `julia +1.10 --project=/tmp/pfw-unitful-110 -e 'using Pkg;
Pkg.test(; julia_args=["--threads=4"])'`. The first isolated-copy attempt lacked
the benchmark fixture included by the existing Lennard-Jones tests; copying that
directory restored the full suite without source changes. `git diff --check`
passes. Root/reference Manifests were preserved.

A separate pinned-reference-environment check loaded the first real ExtXYZ frame,
prepared its design matrix, and confirmed identical training rows from its plain
labels and the same labels with eV/eV/Å units attached (ExtXYZ 0.2.2, Julia 1.13.2;
`/tmp/pfw-unitful-extxyz.log`). Reference inputs were preserved.

## ACE calculator persistence (2026-10-23)

`ACE.ACECalculator(path; kwargs...)` now accepts a resolved TOML export, the
supplied JSON recipe with its hash-checked resolved export, or a native calculator
checkpoint. `ACE.save_calculator(path, calc; parameters=calc.parameters,
state=calc.state)` writes a self-contained, versioned TOML checkpoint. Explicit
parameter/state overrides support training through AtomsCalculators without
installing optimizer results in the mutable calculator. Saving does not mutate
the calculator; one authoritative parameter tree is stored. Complete resolved
basis/reference data, readouts, nonlinear radial/pair scales, numeric types,
array shapes and supported nested state survive loading in a different directory.
Environment, schedule, scratch and derived coefficients are rebuilt with the
requested execution settings. No dependency or core compatibility changes.

The format supports nested NamedTuples/tuples, supported primitive numeric
arrays/scalars, strings and `nothing`; unsupported objects/AD tracers and
nonfinite parameter/state values fail explicitly. A failed validation/write does
not truncate the destination. This is a PotentialFrameWork checkpoint, not an
ACEpotentials JSON export. Recipe import still requires the separately resolved
export; native checkpoints have no reference-file dependency.

All **37 new persistence checks** pass, including changed readouts and nonunit/
zero nonlinear scales, E/F/W and AtomsCalculators force-loss ForwardDiff products
after reloading, Float32 parameters, empty shaped arrays, relocated files, stale
inference caches, installed versus explicit parameter saves, and rejected
shapes/schema versions/nonfinite/tracer/state values. Complete suites pass
**44,463 assertions per process** on Julia **1.13.2 / 32 threads** and
**1.10.12 / 4 threads**, AMD Ryzen 9 9950X. Versions remain AtomsCalculators
0.2.3, SpheriCart 0.2.4, ForwardDiff 1.4.6, Zygote 0.7.14 and TOML 1.0.3.
Logs: `/tmp/pfw-ace-persistence-113-t4.log` (the Julia 1.13 Pkg test subprocess
selected 32 threads despite the parent process's `--threads=4`) and
`/tmp/pfw-ace-persistence-110-t4.log`. Julia 1.10 used the existing isolated
environment `/tmp/pfw-ace-110-flg1rsld` with the current changed source/tests
copied in and its compatible Manifest preserved.

Reproduce with `julia --project -e 'using Pkg;
Pkg.test(; julia_args=["--threads=4"])'`, selecting the desired Julia version
(use the isolated compatible environment for Julia 1.10). A core-only checkpoint
round trip using `angular_backend=:block` also passes without loading SpheriCart
or ACEpotentials (`/tmp/pfw-ace-persistence-optional-free.log`). `git diff --check`
passes. Supplied reference inputs and the existing execution/force ownership
contracts are unchanged.

## CPU ACE phase — current gates (2026-10-22)

M0 reference reproduction passes on Julia 1.13.2 with the supplied, unchanged
reference Manifest and ACEpotentials 0.10.2. All 33 reloaded XYZ structures were
evaluated without fitting. Maximum absolute differences from rounded stored
predictions: energy 2.624960870889481e-8 eV, force 7.195396589665393e-8 eV/Å,
virial 5.642368563485434e-7 eV. Correct per-atom `atom_data.force` provenance is
checked. Full precision recomputed predictions are in
`data/ace_ref/reloaded_predictions.toml`. Reproduce with:
`julia --project=data/ace_ref data/ace_ref/validate_saved_model.jl`.
Saved model SHA256:
`a61a9e57c91cc25663359dd876c47fe9a4c3662ee8d5b4a27a92350c90435f0e`;
XYZ SHA256:
`249c00cc8ce9a2ca6e9d45f74b99c9690a4ade5d00df976deb4fd1365737c62c`.

Actual reconstruction: species order Ti/Al (22/13); spline radial basis with
30 channels, spline pair basis with 12 channels; real spherical L2-normalized
SpheriCart harmonics through l=2; 27 pooled density channels; 167 selected
products; sparse 123×167 invariant map. The 270 supplied parameters are
WB (123×2, column-major) followed by Wpair (12×2, column-major); radial/pair
parameter subtrees are empty in this saved model. Explicit one-body references
are Ti=-1586.0195 and Al=-105.5954 eV. A polynomial recurrence substituted for
these reconstructed splines would change the mathematical oracle.

Fresh pre-edit `Pkg.test()` passes on Julia 1.13.2, 32 test-process threads,
AMD Ryzen 9 9950X (16 cores/32 logical CPUs). Log:
`/tmp/pfw-ace-baseline.log`. Zygote resolved to 0.7.14 in this test environment;
older version/test counts below remain historical. Core Julia compatibility and
the supplied reference Project/Manifest/model/XYZ/generation script are preserved
(their hashes were rechecked at completion). The ignored core Manifest was
resolved for the newly declared standard libraries; its original copy is
`/tmp/pfw-core-manifest-before-ace.toml`. The reference Manifest was never resolved
or overwritten. Julia 1.10 uses an isolated source/environment at
`/tmp/pfw-ace-110-flg1rsld`, not the reference environment.

Implemented gates, in review order:

| Gate | Implementation and evidence |
| --- | --- |
| M1 | `PairPotentials/numerical.jl`: owning numeric inputs, fixed original-order unordered image metadata, explicit E/F/W; legacy callbacks/calculators/SIMD preserved. Edge-local ForwardDiff products and ChainRules cover X/ps E/F/W, including nonlinear LJ. Saved/thunk/zero pullbacks, self-images and deterministic rollout sensitivities pass. |
| M2 | `export_resolved_model.jl` plus `ACE/resolved.jl`: narrow versioned import, source-hash convenience recipe import, all 270 original parameter coordinates, complete spline/transform/envelope and product/coupling tables. 1,136 radial value/derivative assertions pass against the pinned reconstruction. MIT attributions retained under `docs/implementation/licenses`. |
| M3 | `ACE/environments.jl`: full directed flat adjacency and exact image labels, species buckets, safe masked tails, original-order numerical displacements, reusable capacities and new complete write footprints/phases. 28,807 oracle/padding/species/schedule checks pass. Isotopes/custom names are explicitly rejected for element-only channels. |
| M4 | `ACE/serial.jl`: independent scalar radial/angular/density/product/sparse-map/readout reverse, fixed references and imported pair component; AtomsCalculators owning/additive outputs. All 33 fixture E/F/W agree; absolute maxima are 7.275957614183426e-12 eV, 6.163958232718869e-13 eV/Å and 8.446576771348191e-13 eV. Force/strain finite differences, translation/rotation, zero coefficients and empty cases pass. |
| M5 | `ACE/fitting.jl` and optional ACEfit extension: dense original-coordinate E/F/W design, missing-property selection, residual weights, fixed reference subtraction, jointly fitted or fixed pair readout, metadata-derived algebraic prior. ACEfit 0.3.1 weighted rank-deficient prediction recovery and actual QR prior scaling pass (8 pinned-environment checks). The training example's synthetic residual is 3.108624468950438e-14. |
| M6 | `ACE/cpu.jl`: SIMD joint spline Horner evaluation, capacity-managed site scratch, refreshed sparse readout precontraction, exclusive task slots, full-footprint phases and one shared species-sorted force matrix. Actual destinations, widths 1/2/4/8/16, Float32/64 buffers, growth/reset, current-ps updates and multi-site concurrent phases pass. Single-centre phases execute without thread-task creation. No worker-private full force arrays exist. |
| M7 | Optional LuxCore/ChainRules/ForwardDiff extensions: nested explicit readout/radial/pair scales, RNG-based Lux setup/apply, E/F/W rules into positions/parameters (including optional W->X). ACE position and mixed pullbacks share one directional differentiation of analytic reverse stages, without a full Hessian/Jacobian or replay per parameter. Force-loss, saved inputs/workspace mutation, zero/thunk cotangents and pair/ACE rollout X/velocity/parameter sensitivities pass. Periodic multi-image products pass 20 independent finite-difference checks. Full Lux 1.31.4 setup/apply and a nonlinear gradient step pass in `/tmp/pfw-ace-example`. |
| M8 | Optional generated single-input SpheriCart, explicit SIMD angular jets and tiled radial/angular/density fusion are implemented and selectable. All 33 fixtures pass 5,940 CPU width/backend/fusion/thread comparisons; standalone L=0/1/2/4/6 Float32/64 blocks and complete forces pass. SpheriCart 0.2.4's L=0 gradient limitation is handled by the exact constant/zero derivative. Reference backend remains default. Measurement artifacts and reproduction details are in `benchmark/ace`. |

Final package suites pass **44,426 assertions per process** on Julia **1.13.2 /
1 thread**, **1.13.2 / 4 threads**, and **1.10.12 / 4 threads**. Logs are
`/tmp/pfw-ace-final-113-t1.log`, `/tmp/pfw-ace-final-113-t4.log`, and
`/tmp/pfw-ace-final-110-t4.log`. AMD Ryzen 9 9950X, 16 cores/32 logical CPUs.
Core test versions: AtomsBase 0.5.3, AtomsCalculators 0.2.3, SIMD 3.8.1,
SpheriCart 0.2.4, LuxCore 1.5.3, ForwardDiff 1.4.6, ChainRulesCore 1.26.4,
Zygote 0.7.14, Unitful 1.29.0. StaticArrays is 1.9.23 in the Julia 1.13 core
test environment and 1.9.24 in the Julia 1.10/reference environments;
NeighbourLists is 0.6.2 / 0.5.10 respectively. The strict oracle scripts load
ACEpotentials before adding the core project to LOAD_PATH so pinned transitive
versions (including StaticArrays 1.9.24) are already selected.

Optional-free core load/import passes (`/tmp/pfw-ace-optional-free.log`). Reference
scripts are documented in `benchmark/ace/README.md`; numerical logs are
`/tmp/pfw-independent-reference.log`, `/tmp/pfw-all-frame-cpu.log`, and
`/tmp/pfw-ace-periodic-ad.log`. Optimized radial LLVM is saved under
`benchmark/ace/radial-llvm-julia1.13.2.txt` and contains explicit `<4 x double>`
arithmetic; this is not a register-residency claim. `git diff --check` passes.

Final isolated benchmark reports are
`benchmark/ace/results-julia1.13.2-t1.toml` and
`benchmark/ace/results-julia1.13.2-t4.toml`: 191 full-stack records and 72 angular
workload records per process, nine samples after three warmups, BLAS one thread,
no affinity pinning. Controlled synthetic cases vary 8/64 atoms, one/two populated
element channels and spacings 2.8/4.0 Å in addition to fixtures 1/33. Complete
public E/F/W medians (µs) for fixtures 1/33 are 51.920/1144.383 at one thread and
52.060/1151.104 at four threads. The pinned ACEpotentials public calls measure
81.150/1599.675 and 84.390/736.322 µs respectively. Thus there is no general
speedup or multi-thread scaling claim; dense full-site footprints constrain
available independent tasks. The width-4 prepared baseline allocates 64/1728
bytes and stores exactly 48/1296 force bytes for 2/54 atoms, independent of worker
count. Site scratch is 57,448 bytes for one slot and 229,672 for four on these
fixtures. The final force-loss backward median is 126.691 µs with 625,424 bytes
allocated and 1,135,889 bytes of owned saved state on fixture 1 at four threads.
Angular blocks and tiled fusion remain optional: at L=2, 64 neighbours, width 4,
SpheriCart batched/generated-single/explicit-block gradients measure
0.700/1.530/3.340 µs. The complete reports contain all variations and allocation
costs, including slower experiments. `results-prototype-position-vjp-*.toml`
records the superseded repeated-local-gradient implementation only. Reproduction
commands, conventions and limitations are in `benchmark/ace/README.md`.

Remaining limitations: import is deliberately limited to the resolved supplied
component families; nonlinear training exposes diagonal radial/pair scales, not
arbitrary ACEpotentials constructors. Sorting, image selection and cell AD are
not differentiated; numerical position/parameter products use fixed valid
topology. Training targets the explicit primitive/Lux layer, not mutation inside
inference/preparation. Arbitrary higher reverse orders are unverified. Cutoffs
and fixed spline tables are not trainable. Float32 executor buffers use mixed
Float64 imported radial arithmetic. Transform/envelope/interval selection remain
scalar per lane; joint cubic Horner arithmetic is SIMD. Preparation and owning
outputs allocate; there is no Verlet reuse. Dense design construction and saved
AD data allocate. Experimental angular recurrences allocate local temporaries;
tiled fusion retains intermediates for force reversal and need not be faster.
No automatic crossover, general speedup, GPU/MACE, MPI, long-range or NPT/pressure-
coupled differentiable MD support is claimed. All required numerical/derivative
gates listed above pass; these limitations delimit the verified interfaces.

The sections below are historical first-phase records.

## Milestone 0 — decisions

Repository inspected: initially a root module and empty test skeleton. Package
identity, authorship, license and Julia compat `1.10.10` preserved. User handoff
files were untracked on arrival and are retained.

Inspected local dependency sources: AtomsBase 0.5.1 interface/constructors,
SIMD 3.8.1 and 3.7.2 vector operations, NeighbourLists 0.5.10 and 0.6.2
cell_list/types and MIT license (Christoph Ortner, 2018). The 0.6.2 public
`build_cell_list` / `SortedCellList` supplies cell-ID sortperm and CSR ranges,
but requires Julia 1.11; the compatible 0.5.10 uses linked cells. We implement
a focused CPU cell-ID sorting/range scheme independently, reusing that approach
without copying source or depending on undocumented internals. Sparse occupied
cell keys also avoid allocating empty vacuum cells. Spatial sorting is established,
not a research question. No NeighbourLists dependency is needed.

Internal units: Å, eV, eV/Å. AtomsBase positions and cells undergo explicit
Unitful conversion, never blind stripping. Public outputs carry units. Cell
vectors form COLUMNS of H. Only axis-aligned orthorhombic/mixed-periodic and
nonperiodic cells are supported; periodic lengths must strictly exceed twice the
physical cutoff. Isolated infinite diagonal cells are supported on open axes.

Unordered cell pairs owned by lower sorted cell ID; same-cell atom pairs by
lower sorted atom index. Displacement is x_j + H*s - x_i (minimum image).
Callback forces act on central endpoint. Virial is -sum(r_ij ⊗ F_i).
Cutoff is strict r<c and is executor-owned; SmoothPair includes its smooth taper
once in its analytic formula. No skin: every public evaluation rebuilds ordering,
topology and schedule; block geometry is freshly packed during traversal.

Public API: Workspace(; width, precision), evaluate/evaluate!, EnergyOnly and
EnergyForcesVirial requests; cutoff and bind_species extension functions;
PairBlock callback returns per-lane energy and central fx/fy/fz. RadialPotential
adapts U,dU/dr. Production force storage is one shared 3×N sorted matrix. Task
partials hold one energy and nine virial scalars per occupied cell. Original-order
unitful outputs are separate boundary results, not thread-private work arrays.

## Milestones 1–3 — implementation and validation

Implemented `NeighbourList`, `PairPotentials` and `Execution`; root module only
includes/reexports the public API. Scalar W=1 uses scalar floating-point
arithmetic. W=2/4/8/16 uses explicit SIMD.jl vectors. Blocks are consumed
immediately; no global pair-result storage. Inactive lanes have valid indices,
valid species values, safe displacement (1,0,0), and masked outputs. Neighbour
endpoint stores are scalar, handling duplicate padding indices safely. Each
central contribution is retained across its blocks. Model-specific arithmetic
domains still require masks inside custom callbacks. Radial adapters require a
regular origin (or a model-defined singularity policy).

The smooth example supports runtime mixed-species symmetric tables and an
explicit binding stage; no compile-time chemical identities. Independent tests
use a different AtomsBase implementation plus actual FlexibleSystem/isolated
systems, rather than requiring AtomsSystems.

Task ownership is independent of scheduling. Candidate cell adjacency is
sparse, wrapped cells are deduplicated, and a cell-to-colour incidence greedy
algorithm builds phases without a dense conflict matrix. Conservative footprints
include central and all owned candidate neighbour cells. Each cell task stays
serial; tasks in a phase write disjoint cells of the single sorted force matrix.
Barriers separate phases. Per-task partials do not depend on thread IDs. The
schedule version is tied to the freshly rebuilt topology. Direct mutation of
workspace internals is not a supported invalidation API.

Checks include independent O(N²) energy/force/virial agreement; finite differences
of forces and homogeneous dilation/shear; total force, translation, permutation;
Float32/Float64 and widths 1/2/4/8/16; empty/single/isolated systems; cutoff edges,
coincident smooth interactions and inactive lanes; explicit unit conversion;
unknown species errors; periodic faces, multi-box offsets and boundary crossings;
small wrapped grids with lengths 4.01/6.1/8.1/11 Å at cutoff 2 Å; mixed/open axes;
changed positions/species/cell and repeated output resets. Tests enumerate every
active unordered pair and assert exactly one owner against the oracle. They
check phase disjointness and actual force destinations against footprints, not
just central-cell ownership. Unsupported cells/cutoffs and stale schedule
versions are rejected.

Generated-code inspection: `benchmark/inspect_simd.jl` on Julia 1.10.12 produced
LLVM `fmul <4 x double>`, `fadd <4 x double>` and `fsub <4 x double>` in the W=4
callback. W=1 is a genuine scalar baseline. This establishes vector arithmetic;
it is not evidence about whole-driver native spills or a GPU execution path.

## Milestone 4 — reuse and measurement method

Coordinate/species/permutation vectors and force/partial arrays are reused when
sizes permit. Public calls rebuild temporary extraction, cell keys, adjacency,
footprints and schedules, and rebind species every time. There is no cached
Verlet list. Geometry is freshly packed during the fused evaluation/accumulation
loop. A fixed-topology numeric entrypoint is available internally for benchmarks,
with unchanged system/parameters as its precondition; it resets all outputs.

Allocation inspection found a boxed self-captured virial tuple in the first
implementation. Moving its update into an inline function removed the box.
A warmed fixed-topology serial 2048-atom W=4 diagnostic then allocated zero
bytes (measured with `@allocated`). Benchmark wrappers can include a small dynamic
call overhead; the recorded raw allocation figures are authoritative for those
calls. Parallel execution allocates Julia scheduling tasks/barriers per phase;
these are not atom-by-thread force arrays.

`benchmark/run.jl` measures seven samples after three full warmups, reporting
minimum/median/maximum complete times and stages, Julia allocated bytes, workspace
size, force/partial storage, phases/tasks/work per phase, task imbalance and
lane occupancy. Fresh geometry and callback/accumulation are fused; their combined
stage is timed, plus a separate diagnostic geometry/occupancy scan. That scan is
not an additive estimate of callback time. Partial merge and output permutation
are timed separately. Fixed-topology timings amortise preparation over all samples;
public timings include every rebuild. No per-thread-force baseline, RSS/peak
memory measurement, native spill analysis or pinned-affinity study is claimed.

## Milestone 5 — handoff

README documents input accessors, units, force/virial conventions, scalar/SIMD
callback contract, species binding, radial adapter, widths/precision, threaded
usage, ownership, scheduling, lifetime and geometry restrictions. The runnable
custom callback example includes full Cartesian and radial interfaces and an
energy-only request. Benchmark scripts and reproduction instructions are included.
Existing CI Julia versions/platforms are retained, adding 1/4-thread jobs.

Final validation is recorded in [TEST-RESULTS.md](TEST-RESULTS.md): **1,790 passing
assertions per process** on Julia 1.10.12 with 1, 4 and 32 threads, and Julia
1.13.2 with 4 threads. The four-thread Julia 1.10 run used `Pkg.test()` in its
isolated test environment. The exact 1.10.10 patch and prerelease Julia have not
been run locally. Examples and representative SIMD LLVM inspection passed.

AtomsCalculators 0.2.3 public interface was inspected locally. Its vector of
static-vector forces differs from this package's 3×N boundary result. A small
boundary adapter is follow-up work rather than an unverified integration claim.
No ACE/MACE, GPU, long-range, training, MPI, per-atom virial/stress or NUMA
few-buffer implementation is included. Extension boundaries remain: site-level
complete neighbourhood preparation for ACE; layer-wide feature stages/reverse
accumulation for MACE; separate GPU mapping/synchronisation and long-range models.

All gates through milestone 5 have implementation and correctness evidence.
Performance limitations are recorded below; there is no speedup acceptance gate.

### Recorded benchmarks — 2026-10-20

AMD Ryzen 9 9950X, 16 cores / 32 logical CPUs, Linux x86_64, generic Julia
compiler target. Julia 1.10.12; AtomsBase 0.5.3, SIMD 3.8.1, Unitful 1.29.0;
Float64. Separate 1/4-thread processes; no CPU affinity control. Three complete
warmups and seven measured samples per case; no concurrent tests/benchmarks.
[Full summary](../../benchmark/SUMMARY.md) links raw records for both thread counts,
including uniform/clustered, single/mixed species, widths 1/4/8 and 128/2048/8192
atoms (8192 only in the four-thread process).

Selected mixed-species complete-call medians in the four-thread process
(serial modes use one worker), milliseconds:

| N / distribution | Scalar serial | W=4 serial | W=4 phased parallel |
|---|---:|---:|---:|
| 128 uniform | 0.087 | 0.087 | 0.135 |
| 2048 uniform | 1.733 | 1.748 | 1.519 |
| 2048 clustered | 6.027 | 4.555 | 4.826 |
| 8192 uniform | 7.181 | 7.326 | 6.291 |

For uniform N=8192, W=4 phased execution: public min/median/max
6.134/6.291/8.918 ms; allocated 9,253,088 bytes/public call. Median stages (ms):
extraction 0.472, sorting/SoA packing 0.623, candidate topology 3.390, schedule
0.582, parameter binding 0.00026, fused fresh geometry/callback/accumulation
1.154, partial merge 0.011, final permutation 0.037. Stage medians do not
necessarily sum to the median complete call. A separate serial geometry-only
scan including occupancy counting took 1.173 ms median. The warm fixed-topology
parallel traversal took 0.720 ms median and allocated 69,984 bytes; it benefits
from prebuilt/hot metadata and is not the complete public cost.

Force storage was exactly 3×N×8 bytes: 3,072 / 49,152 / 196,608 bytes for
128/2048/8192 atoms, independent of worker count. N=2048 force storage was
49,152 bytes in both 1- and 4-thread processes. Uniform N=8192 total retained
workspace was 2,053,832 bytes, of which energy/virial task partials were 279,840
bytes. Boundary output is one additional original-order 3×N result. Public
construction allocations remain much larger than retained workspace memory.

Uniform N=8192 used 31 phases; occupancy ranged from 190 tasks down to one task
per phase. Only 8.5% of W=4 candidate lanes were physically active (4.5% for
W=8); cell tails and the conservative candidate stencil account for wasted work.
Clustered N=2048 had 24 phases but task pair work ranged 0–52,031, exposing
severe imbalance. Raw records include complete per-phase work/occupancy vectors.

Limitations observed: scalar is competitive for sparse uniform systems; SIMD
helps the dense cluster but wider vectors are not universally better. Small
systems lose time to phase barriers, and dense clusters can remain dominated by
one serial cell task. Serial topology construction dominates larger complete
calls, limiting benefits from faster threaded traversal. Next performance work
should target topology allocation/reuse, candidate packing and phase utilisation;
cache-line false sharing is possible and not yet measured. No universal speedup,
thread-private-force comparison, peak RSS, affinity tuning or native spill result
is claimed. These regressions are recorded rather than hidden by kernel-only
numbers or weakened scheduling checks.


## Follow-up — ChemicalSpecies representation

At the user's request, species storage and callback fields now use AtomsBase
`ChemicalSpecies` rather than compact integer species IDs. This supersedes the
initial specification's compact-ID representation. `Workspace.species` and
`species_map` are concrete `Vector{ChemicalSpecies}`; blocks carry one central
ChemicalSpecies and an NTuple of neighbour ChemicalSpecies values. System access
still occurs only at extraction; no chemical species are encoded in types.

`SmoothPair` normalises symbol labels to ChemicalSpecies and owns a concrete
`Dict{ChemicalSpecies,Int}` mapping to its private parameter-table indices. These
indices describe table rows, not atom species. Binding happens once per call;
the read-only dictionary is shared safely by cell tasks. Species identity uses
strict `isequal`/`hash`, preserving isotope and atom-name fields. AtomsBase `==`
is a wildcard matching relation and is unsuitable for this lookup. The minimum
AtomsBase version is now 0.5.3, which supplies the strict equality/hashing API;
Julia compatibility remains 1.10.10.

The independent oracle and fixtures now use ChemicalSpecies and strict identity.
New tests cover species field types, isotope/name preservation, distinct table
entries, canonical duplicate labels, missing exact species parameters, permutation,
and radial callbacks typed on ChemicalSpecies, for scalar and SIMD widths.
Earlier benchmark records above and in benchmark/SUMMARY.md describe the previous
integer-species implementation and are retained as historical measurements.

Validation after this change: **1,837 assertions pass per process**, including 47
new ChemicalSpecies assertions, on Julia 1.10.12 with 1 and 4 threads (Pkg.test)
and Julia 1.13.2 with 4 threads. Julia 1.10 tests used a temporary development
environment to avoid changing the existing Julia 1.13 Manifest.toml. The callback
example and updated SIMD LLVM inspection passed; W=4 still emits vector floating
point operations. The full historical benchmark matrix has not been rerun for
this representation; dictionary lookup costs should be measured before comparing
its performance against the previous integer-species implementation.
A warmed serial W=4 fixed-topology traversal of 2,048 mixed-species atoms still
allocated **0 bytes** with the ChemicalSpecies lookup (Julia 1.10.12, four-thread
process, measured using `@allocated` inside a specialised function).


## Follow-up — AtomsCalculators interface

Pair potentials now implement the AtomsCalculators 0.2.3 interface directly.
Inspected the installed interface, interface-generation helpers, low/high-level
example and conformance tests, plus the upstream repository documentation:
https://github.com/JuliaMolSim/AtomsCalculators.jl . AtomsCalculators 0.2.3 supports
Julia 1.6+ and AtomsBase 0.5, compatible with this project's Julia 1.10.10 minimum.
AtomsCalculators and StaticArrays are direct dependencies.

`AbstractPairPotential` is the opt-in extension point for custom concrete
callbacks; SmoothPair and RadialPotential subtype it. `CalculatorInterface`
implements potential_energy, forces/forces!, virial, energy_forces/energy_forces!,
energy_forces_virial/energy_forces_virial!, unit/zero/type helpers and individual
or tuple calculate requests. Combined requests run one traversal (tested through
workspace topology versions). The low-level matrix-output evaluate API remains
available. AtomsCalculators gets Vector{SVector{3,...}} forces and SMatrix virial,
with eV / eV per Å units and the same original atom order and virial convention.

The optional `neighbourlist` keyword accepts the existing NeighbourList.Workspace;
otherwise a workspace is constructed inside the call. A supplied workspace is
refreshed/rebuilt using the current system and cutoff, never treated as a stale
cached topology. Width and precision arguments must match a supplied workspace.
Automatic precision follows Float32 input coordinates or defaults to Float64;
empty systems default to Float64 unless overridden. Default zero/type utilities
also work for empty systems. Generic system access remains outside hot loops.

Important upstream contract: forces! and both combined ! methods are ADDITIVE,
not replacement writes. The existing sorted workspace is cleared each evaluation,
then its original-order increment is added into the caller's destination. Length
mismatches raise DimensionMismatch. Boundary conversion currently allocates an
additional vector of static vectors; mutating methods are not advertised as
allocation-free. No per-thread force arrays, atomics or new pair-result storage
were introduced. Previous benchmark results do not include this boundary adapter.

The calculate interface accepts positional or keyword parameters/state. Potential
parameters are bundled; nonempty parameter overrides are explicitly rejected,
while nothing/empty NamedTuple select bundled values. State passes through
unchanged. Extra keywords are accepted for AtomsCalculators composition, as
required by its conformance tests. No training/stateful model API is claimed.
This completed adapter supersedes the earlier AtomsCalculators deferral above.

Validation: **2,332 assertions pass per process** under Pkg.test on Julia 1.10.12
with 1 and 4 threads, and Julia 1.13.2 with 4 threads. AtomsCalculators 0.2.3,
StaticArrays 1.9.23, AtomsBase 0.5.3, SIMD 3.8.1 and Unitful 1.29.0 were used.
The 495 added assertions cover upstream conformance for SmoothPair and radial
adapters (118), precision/unit/output shape and workspace lifetime (336), and
additive destinations, individual/tuple calculate, state/parameter handling and
energy-only callbacks (41). Core numerical, periodic ownership and parallel write
footprint checks remain enabled. Tests compare component values explicitly:
Unitful vectors nested in static vectors require component-wise or stripped-unit
comparisons for tolerances, rather than Julia's generic nested-array norm.

The 1.10 test environment was resolved separately with a compatible PrecompileTools;
the ignored root Manifest was resolved under Julia 1.13 to add AtomsCalculators
and its dependencies, without upgrading existing package versions. The updated
custom callback example passed on Julia 1.13.2 with four threads. `git diff
--check` passed. No speedup or new benchmark timing is claimed for this adapter.

## Follow-up — Lennard-Jones comparison setup

Inspected EmpiricalPotentials 0.2.4 and 0.2.6 LennardJones,
AtomsCalculatorsUtilities 0.1.8/0.1.9 site/pair assembly, NeighbourLists
0.5.10/0.6.2 and AtomsSystems 0.1.2 public constructors.
The reference LJ is shifted/tilted in s=rmin/r, not in r. The benchmark-only
ArgonLJ callback matches U=epsilon*(f(s)-f(sc)-(s-sc)*f'(sc)), f=s^12-2s^6,
and uses its analytic Cartesian derivative; upstream uses ForwardDiff. This
matches the same model without modifying the reference's implementation.

Argon example values: epsilon=0.010323 eV, sigma=3.405 Å,
rmin=2^(1/6)*sigma, cutoff=2.5*sigma, reference emin=-epsilon. Deterministically
perturbed FCC systems use a=5.26 Å, 0.04 Å component perturbations and
N=4*m^3. AtomsSystems.generic_system constructs the systems; periodic positions
are wrapped into their original cell. Periodic m>=4 satisfies the framework's
minimum-image cutoff restriction. Parameters are illustrative, not a new fit.

AtomsSystems and EmpiricalPotentials are test-only extras, not runtime package
requirements. The benchmark has a separate environment pinning the inspected
reference versions and a relative-path framework dependency. It launches one
fresh Julia process per thread count, sequentially to avoid cross-case contention.
Complete public energy/forces/virial calls rebuild lists every time in both
packages; workspace buffers are reused where supported. Separately labelled
fixed-topology diagnostics exclude setup and require unchanged geometry. They
are not public Verlet reuse or a complete-call speedup claim.

Reference serial mode uses Folds.SequentialEx and ntasks=1; threaded mode uses
Folds.ThreadedEx and ntasks=Threads.nthreads(). The recorded benchmark targets
EmpiricalPotentials 0.2.6 / NeighbourLists 0.6.2 sort-based construction on Julia
1.13.2, keeping its default CPU backend. Julia 1.10 can resolve the older 0.5
backend, whose list construction is serial. Site assembly has full per-task force partials;
this is the external baseline, not a new production accumulation mode here.
The framework uses W=1 scalar, W=4/8 SIMD, serial or phase-parallel execution,
with one shared sorted force matrix. Both paths return the same AtomsCalculators
boundary units/order/shapes. Every benchmark case checks energy, forces and
virial agreement before and after timing; failure aborts the benchmark.

### Lennard-Jones checks and recorded measurements

The complete suite passes **2,616 assertions** per process, including 284 new LJ
checks, on Julia 1.10.12 with 1 and 4 threads and Julia 1.13.2 with 4 threads.
All use EmpiricalPotentials 0.2.6; Julia 1.10 resolves
AtomsCalculatorsUtilities 0.1.8 / NeighbourLists 0.5.10, while Julia 1.13 uses
0.1.9 / 0.6.2. Regression cases include two-atom distances at/below/above cutoff,
open FCC clusters, mixed/fully periodic FCC, small wrapped grids, widths 1/4/8,
serial/threaded agreement, force conservation, caller immutability and a separate
energy finite difference. Real coincident LJ atoms raise DomainError.

The recorded performance matrix uses **Julia 1.13.2**, EmpiricalPotentials 0.2.6,
AtomsCalculatorsUtilities 0.1.9, NeighbourLists 0.6.2, AtomsSystems 0.1.2,
AtomsBase 0.5.3, AtomsCalculators 0.2.3, ForwardDiff 1.4.6, Folds 0.2.10,
SIMD 3.8.1, StaticArrays 1.9.23 and Unitful 1.29.0. Hardware: AMD Ryzen 9 9950X,
16 cores/32 logical CPUs; compiler target znver5, Float64, uncontrolled affinity.
Four independent, sequential worker processes use 1/2/4/8 Julia threads; sizes
256/864/4000/16384 atoms; 3 warmups and 7 measured samples. No concurrent tests
or benchmarks ran during this matrix. Both full rebuild and fixed-topology
regimes cover reference serial/threaded modes, W=1 scalar and W=4/8 serial/phased
modes: **224 measured cases**, all passing pre/post numerical checks.

Largest recorded absolute deviations across all cases: energy 5.912e-12 eV,
force component 4.025e-16 eV/Å, virial component 5.059e-12 eV. The matching cutoff
formula is important; an ordinary distance-shifted LJ would not pass this check.

Selected complete-call medians (reference threaded site executor vs framework
W=4 phased), milliseconds; all include fresh topology construction:

| N | Threads | EmpiricalPotentials | PotentialFrameWork | EP/PFW median ratio |
|---|---:|---:|---:|---:|
| 256 | 1 | 2.798 | 0.221 | 12.64 |
| 864 | 8 | 7.013 | 1.516 | 4.63 |
| 4000 | 8 | 10.445 | 1.551 | 6.74 |
| 16384 | 1 | 171.789 | 22.523 | 7.63 |
| 16384 | 8 | 38.132 | 5.749 | 6.63 |

For N=16384, 8 threads, W=4: complete min/median/max 5.711/5.749/6.255 ms,
versus reference 27.179/38.132/63.380 ms. Median framework stages (ms): extraction
0.153, sorting/packing 0.741, topology 0.470, schedule 0.086, binding 0.001,
fused geometry/callback/accumulation 4.243, partial merge 0.002, output permutation
0.033. Static-vector boundary conversion and timer overhead are additional;
separate stage medians need not sum to the median call. The fixed-topology
medians are 4.832 vs 15.823 ms, excluding setup and assuming unchanged geometry.

Complete-call allocated bytes at that size/count: framework 3,304,544 vs reference
186,042,456 (Julia GC-visible bytes, not peak memory). The framework sorted force
workspace is exactly 393,216 bytes at every thread count; task scheduling adds
small thread-dependent allocations. Reference allocations include per-site
neighbour/gradient arrays and full per-task force partials/reduction; the allocation
gap is not attributed solely to its thread-private force storage.

Measured limits: m=4 and m=6 produce 2³/3³ occupied cells, with one task in each
of 8/27 phases, so the framework cannot parallelise their accumulation. Cached
reference evaluation at N=864, 8 threads is **0.594 ms vs our 1.470 ms**; the
reference wins there. At N=16384 the framework has 27 phases, initially 36 tasks
per phase and a tail down to one task. Complete-call 1-to-8-thread scaling is
3.92× here vs reference 4.51×; fixed-topology scaling is 4.35× vs 5.22×. Wider
vectors are not a win in serial: at N=16384, W=1/4/8 complete serial medians are
22.284/22.186/22.524 ms. These data show no material SIMD gain for this workload;
they do not isolate its cause. Geometry packing/scatters/reduction are candidates
for subsequent profiling, not assumed bottlenecks.

These ratios concern whole implementations: a single-species analytic callback
with precomputed cutoff constants versus a general multi-species LJ evaluator,
AD derivatives, directed vs unordered pair ownership, topology construction,
allocations and scheduling. They are not a SIMD-only or universal speedup claim.
No production optimisation or weaker correctness check was introduced to improve
the benchmark. CPU affinity/NUMA, peak RSS, native spill and cross-Julia-version
performance studies remain unmeasured.

Reproduction, complete tables/variation, source/manifest, raw CSVs, exact package
metadata, phase occupancy log and PNG/SVG scaling plots are in
[benchmark/lennard_jones](../../benchmark/lennard_jones/README.md) and
[RESULTS.md](../../benchmark/lennard_jones/RESULTS.md). Reports were regenerated
and the plot visually inspected; `git diff --check` passes.

## Follow-up — explicit W=8 scaling and Zen 5 code inspection

The recorded LJ matrix already includes W=8 at every size/thread count, both
serial and phased, rebuild and fixed-topology (64 W=8 cases with pre/post numerical
agreement checks). No new timing samples were substituted. Reporting now shows
W=8 ranges, W4/W8 ratios, W1/W4/W8 serial ratios and independent W=8 thread
scaling. Both widths appear in the main time plots, and a separate ratio plot
shows W=4 time / W=8 time against thread count for each size. Ratios >1 favour
W=8; they are ratios of measured medians, not confidence intervals.

Added `benchmark/lennard_jones/inspect_widths.jl`, inspecting both the callable LJ
kernel and the actual fused cell-task traversal at W=1/4/8. Executed on Ryzen
9950X with Julia 1.13.2, `-C native`, target znver5, Float64. Native code confirms
W=4 callback packed-double arithmetic uses YMM registers and W=8 uses ZMM, e.g.
`vmulpd zmm3, zmm0, zmm0`. W=8 processes eight Float64 lanes (512 bits), and its
fused traversal also contains ZMM arithmetic. LLVM/native output and a generated
inspection summary are retained in `benchmark/lennard_jones/codegen`.

Static packed-double arithmetic-line counts: W=4 callback 32 YMM / 0 ZMM;
W=8 callback 0 YMM / 32 ZMM; W=8 traversal 0 YMM / 60 ZMM. These are source-code
inspection counts, not executed instruction counts or throughput. Compiler
vectorisation also emits some ZMM arithmetic in W=1/W=4 traversal, so a logical
pair width of one does not mean every machine operation must be scalar.

AVX-512 is in use, but this complete workload has no broad W=8 speedup. At
N=16384, W=4/W=8 phased rebuild medians by threads are:
1: 22.523/22.658 ms; 2: 13.538/13.421; 4: 8.462/8.841; 8: 5.749/6.219.
Eight-thread fixed-topology medians are 4.832/4.941 ms. W=8 complete-call scaling
from one to eight threads is 3.64×, versus 3.92× for W=4. At one thread, serial
W=1/4/8 medians are 22.284/22.186/22.524 ms. The data do not establish why wider
vectors fail to improve the complete path; geometry packing, masks/lane occupancy,
endpoint updates, tensor reduction and scheduling remain profiling candidates.
No default width, physics, accumulation or scheduling code was changed.

Validation for this reporting/code-inspection change: regenerated tables and
both PNG/SVG plots from all 224 existing CSV cases; inspected the ratio plot and
updated time plot visually; successfully generated specialised callback and
traversal LLVM/native code at all three widths; `git diff --check` passes. The
previous full suite and benchmark numerical checks remain the evidence for the
unchanged evaluation implementation; the full timing matrix was not rerun.

## Energy derivatives with respect to potential parameters (2026-10-20)

The requested scope now includes energy parameter AD through AtomsCalculators.
This supersedes the earlier bundled-only parameter restriction. Explicit nonempty
parameters are accepted by `calculate(Energy(), sys, p, parameters[, state])` and
its keyword/energy-only tuple forms. The result retains unitful energy and the
unchanged caller state. Explicit force/virial requests fail clearly; they do not
silently drop derivative information. Passing explicit parameters to other public
property helpers also errors, directing callers to `calculate(Energy(), ...)`.

Parameter contract: a flat NamedTuple of real scalar/vector fields;
`AtomsCalculators.get_parameters(p)` exposes values and `with_parameters(p, ps)`
reconstructs a callback without mutating the template. The physical cutoff must
remain fixed and is checked against the template. System geometry, species,
sorting, topology and ownership are outside this derivative. SmoothPair exposes
one independent amplitude per symmetric upper-triangle entry, in column order;
its off-diagonal parameter changes both entries. Its amplitude and radius types
are now independent, so changing parameter precision does not round the fixed
cutoff. `set_parameters!` follows AtomsCalculators' permitted return-new-object
contract. Generic custom models can extend the same hooks; the benchmark ArgonLJ
now exposes epsilon/rmin, retaining its independent fixed physical cutoff and
recomputing both inverse-distance cutoff coefficients with parameter arithmetic.

Implementation: ParameterEnergy is a submodule. Numeric geometry/topology use the
existing NeighbourList Workspace and Execution.geometry ownership/masks. Only
active lanes invoke scalar, parameter-promoted energy callbacks; Dual values never
enter SIMD.Vec or floating-point workspace output storage. Task-local scalar
energy partials (possibly Dual) are separate from force storage. Serial and phased
parallel evaluation use the same candidate traversal; no full force array per
thread is introduced. All public parameter evaluations rebuild geometry, including
calls made by a gradient backend. The ordinary SIMD energy/force/virial executor
is retained. Callback parameter arithmetic must be generic over scalar real types.

ForwardDiff and ChainRulesCore are optional weak dependencies with package
extensions (compat ForwardDiff 0.10.36/1, ChainRulesCore 1.25+). Direct ForwardDiff
energy differentiation preserves dual numbers, including nested Dual Hessians.
The ChainRulesCore extension defines an rrule on numeric-eV `parameter_energy`;
`calculate` wraps its result in eV outside the primitive. This avoids differentiating
mutable sorting/workspace construction. System and bundled-model cotangents are
explicitly NotImplemented, rather than false mathematical zeros; differentiate
the explicit parameter argument. The `energy_parameter_gradient` hook allows
analytic local gradients. SmoothPair uses exact unit-basis energies for its linear
amplitudes without ForwardDiff; the ForwardDiff extension supplies a generic
fallback for other models. Thus reverse mode for a generic potential presently
uses local forward AD. Gradients are captured when constructing the pullback, so
subsequent changes to a workspace, system or parameter vector cannot corrupt it.
Zero/thunk cotangents and parameter projection are covered.

Validation includes independent all-pairs mixed-species amplitude derivatives;
W=1/4/8 and serial/phased execution; scalar/vector parameter fields; Float32
geometry and parameter precision; empty/single-atom and no-active-pair cases;
fixed-cutoff and invalid-parameter errors; direct rule scaling/zero/thunk inputs;
repeated pullbacks after input mutation; and actual Zygote gradients through
positional and keyword AtomsCalculators.calculate. LJ tests compare to analytic
epsilon/rmin derivatives including precomputed-cutoff derivatives, finite
differences on periodic FCC Argon, and direct ForwardDiff/Zygote agreement.
ForwardDiff Hessians are also checked. The runnable example gives E=0.0421875 eV
and the same amplitude gradient [0.421875] from ForwardDiff and Zygote. A separate
process tested the exact SmoothPair rrule with ChainRulesCore loaded and the
ForwardDiff extension absent; it passed.

Limitations: first-order reverse energy parameter gradients only; no reverse
Hessian, system/cutoff derivatives, force/virial parameter AD or training loop.
Scalar parameter callbacks do not inherit the earlier W=4/W=8 SIMD performance
measurements. Generic local ForwardDiff scales with parameter count; the exact
SmoothPair reference gradient currently does one rebuild/evaluation per amplitude.
Neither approach is claimed suitable for large parameter models. Custom analytic
or backend-specific gradient implementations can replace this hook. No AD speedup
or timing measurement is claimed and the full LJ performance matrix was not rerun;
its existing measurements describe the ordinary, unchanged evaluation path.

Final regression results: **2704 assertions passed** on Julia 1.13.2 with one
thread and Julia 1.10.12 with four threads (full package suites, including existing
numerical, ownership, write-footprint and AtomsCalculators conformance checks).
An earlier Julia 1.13.2 four-thread run also passed before the final cutoff-type
precision checks were added. Final suites used AtomsBase 0.5.3, AtomsCalculators
0.2.3, ForwardDiff 1.4.6, ChainRulesCore 1.26.4, Zygote 0.7.13,
EmpiricalPotentials 0.2.6 and SIMD 3.8.1 on the Ryzen 9 9950X. Julia compatibility
remains 1.10.10; tests used installed 1.10.12. Logs are
`/tmp/pfw-ad-tests-1-final.log` and `/tmp/pfw-ad-tests-1.10-final.log`.
The final runnable example and standalone ChainRules-only check passed;
`git diff --check` passed. No unrelated changes were reverted.

## General cells, multiple periodic images and small-system traversal (2026-10-21)

The requested scope now includes triclinic/rotated cells and all images within
cutoff. This supersedes the earlier orthorhombic/minimum-image restriction.
The public AtomsBase boundary and column-vector cell convention are unchanged.
The axis-aligned orthorhombic minimum-image path remains for periodic lengths
strictly above twice the cutoff; otherwise images are enumerated explicitly.

General-cell construction sorts fractional-coordinate cell keys into sparse
occupied ranges. If B is the internal basis, reciprocal row norms give the
conservative bounds `|delta_fractional[a]| < cutoff*norm(inv(B)[a,:])`.
Periodic bin counts use the corresponding cell heights. Stencil radii follow
these bounds, including when one cell contains several images of an atom.
Wrapping a bin carries an explicit integer image shift; different shifts of the
same occupied cell are retained. Geometry is `x_j + B*s - x_i`, packed directly
into existing SIMD blocks. Fractional component rounding is never assumed to
find the nearest image of a skewed lattice. No dependency source was copied.

Finite independent periodic generators are required. For mixed periodicity,
open vectors are replaced internally by a Cartesian completion of the periodic
subspace; their supplied zero/infinite lengths need not define an invertible
full cell. Original periodic generators and physical Cartesian positions are
preserved modulo lattice translations. Fully open systems retain Cartesian bins.

Ownership: sorted cell pairs are ordered as before; for distinct atoms, all
image shifts with the lower sorted atom first are considered. For self images,
only lexicographically positive shifts are owned. The zero-shift self pair is
excluded. This equals half the directed periodic sum: self-image energies and
virials contribute once per ±image pair and endpoint forces cancel. Distinct
images can update the same atom; scalar lane endpoint updates remain safe.
Image entries and unique cell adjacency are separate: scheduling deduplicates
physical destination cells, while evaluation retains distinct image shifts.
Production still uses one shared sorted force matrix, with scalar/tensor task
partials and validated disjoint write footprints per phase.

`Workspace(bruteforce_threshold=32)` uses a single all-atom range at or below the
threshold. `evaluate` and AtomsCalculators accept the same keyword; a supplied
workspace owns its setting and conflicting keywords error. Zero forces cell
traversal for nonempty systems. The brute path avoids spatial neighbour search
and uses one serial task, even when parallel execution is requested. It shares
SIMD evaluation, physical cutoff, image ownership, accumulation and parameter
AD traversal with the cell path. For wrapped fractional coordinates, image
indices satisfy `|s[a]| < reach[a]+1`, giving integer bounds `ceil(reach[a])`.
The default 32 is provisional, not a measured universal performance crossover.

Dependency investigation: installed AtomsSystems 0.1.2 `src/utils.jl`
`distance_vector`/`distance` accepts atom ranges and system pairs, but returns one
component-wrapped displacement per pair. It does not enumerate multiple images;
component wrapping also cannot establish a nearest-image result for arbitrary
skewed bases. An independent internal implementation supplies the required
image topology, SIMD layout and spatial task footprints without making
AtomsSystems a required concrete system or runtime dependency. The previously
inspected NeighbourLists sorting/range approach continues to underlie storage.

Validation adds an independent original-order, exhaustive fixed-lattice-cube
oracle (shell convergence checked separately), not the production reciprocal
stencil. Coverage includes skewed cells where fractional rounding fails, rotated
and left-handed generators, fully/mixed periodic cells, tiny cells with several
images and self images, exact self-image cutoff, unwrapped translations,
permutation invariance, empty systems, Float32/Float64, widths 1/4/8, brute/cell
agreement, actual pair-image ownership and actual force-write footprints.
Force and shear/dilation virial finite differences include multi-image cells.
ForwardDiff and Zygote parameter derivatives include self-image contributions.
Repeated phase-parallel execution on a 257-atom skewed system checks phases with
multiple concurrent tasks; workspace reuse switches geometry and traversal paths.

Final Julia 1.13.2 four-thread Pkg.test passes **4,767 assertions** on AMD Ryzen 9
9950X (16 cores/32 logical CPUs). Logs: `/tmp/pfw-final-113-t4.log`.
AtomsBase 0.5.3, AtomsCalculators 0.2.3, AtomsSystems 0.1.2, SIMD 3.8.1,
StaticArrays 1.9.23, Unitful 1.29.0, EmpiricalPotentials 0.2.6,
ForwardDiff 1.4.6, ChainRulesCore 1.26.4 and Zygote 0.7.13 were used.

Limitations: no skin/Verlet reuse, lattice reduction or performance crossover
calibration. Highly skewed or very small periodic cells may generate many
candidates and conservative phases may limit concurrency. No all-pair distance
matrix or per-thread force matrix is allocated. Existing LJ benchmark timings
predate this change and are not performance evidence for general-cell/image or
brute-force traversal. No new speedup is claimed. Package identity, dependency
requirements and Julia compatibility remain unchanged.

The final source also passes **4,767 assertions** on Julia **1.10.12 with four
threads**, resolving dependencies in an isolated `/tmp/pfw-110` source snapshot
so the root Julia 1.13 Manifest is preserved. Compatible PrecompileTools was
downloaded into the temporary depot; no package requirements were changed.
Direct/test dependency versions match those listed above; Julia 1.10 resolves
NeighbourLists 0.5.10 / AtomsCalculatorsUtilities 0.1.8 for the external LJ
reference. Log: `/tmp/pfw-final-110-t4-rerun.log`. An earlier one-thread Julia
1.13.2 run passed the same 4,767 checks before tightening the conservative
brute-image bound. Final `git diff --check` passes.

## User neighbour-strategy benchmark (2026-10-21)

Added exported `benchmark_neighbourlist(sys, cutoff; width=4,
precision=Float64, parallel=false, samples=20, warmup=3, io=stdout)` in the
NeighbourBenchmarks submodule. Numeric cutoffs mean Å; unitful lengths convert
explicitly. The helper creates two private workspaces, forces each strategy,
warms the same timed wrapper and reusable buffers, then alternates measurement
order. It does not change the input system, existing workspaces or defaults.

Each sample measures complete preparation/scheduling plus SIMD candidate geometry
traversal and active pair-image counting, using Execution.geometry and the same
image ownership/ranges/phases as evaluation. Threaded scans use task-local integer
counts and the existing phase barriers; brute force retains its serial task.
Both methods must produce identical active pair-image counts throughout the run.
There is no all-pair distance matrix or per-thread force matrix. No new dependency
or dependency API is introduced; Base.@timed supplies time/allocation samples.

The printed report gives min/median/max total times, median preparation/traversal
times, median allocated Julia-visible bytes, active image count, settings, Julia
version and hardware. The returned NamedTuple preserves raw seconds/bytes samples,
method summaries, the lower median geometry method and cell/brute median ratio.
Empty systems report equivalent strategies. Output can be suppressed with
`io=nothing`. The report shows how to configure either Workspace threshold.

This is explicitly a geometry comparison: system + cutoff cannot measure an
unspecified potential. Callback evaluation, force/virial accumulation and output
conversion are excluded. Real-potential complete-call timing is needed when
those costs influence the choice. GC remains enabled; interleaved samples reduce
order bias but do not eliminate system noise. Stage medians need not add up to
the total median. No global threshold calibration or automatic setting mutation
is performed. README and benchmark documentation include usage and these limits.

New checks compare reported image counts with the independent exhaustive oracle
for empty/single-atom, open, periodic and skew multi-image systems. They verify
unchanged inputs, unit conversion, Float32/Float64 and widths 1/4, serial/phased
scans, returned samples/summaries, printed scope, silent operation and invalid
sample/warmup/cutoff options. Timing ratios are never asserted as performance
requirements.

Final validation: **5,010 assertions passed** per full Pkg.test suite on Julia
1.13.2 / 1 thread and Julia 1.10.12 / 4 threads (243 added benchmark checks).
Logs: `/tmp/pfw-choose-113-t1.log`, `/tmp/pfw-choose-110-t4.log`. The Julia 1.10
run uses the isolated source snapshot/environment described above. Direct/test
dependency versions are unchanged from the general-cell validation. Final
`git diff --check` passes.

A standalone run after both suites finished used Julia 1.13.2, AMD Ryzen 9 9950X,
Float64, width 4, 1 thread, serial execution, uncontrolled affinity, 20 samples
and 3 warmups. Reproduction from the repository root:

```sh
julia --project --startup-file=no --threads=1 -e 'using PotentialFrameWork, AtomsBase, Unitful; include("test/fixtures.jl"); for n in (32,256); benchmark_neighbourlist(fixture(n;L=12.),2.;samples=20,warmup=3); end'
```

Observed total geometry min/median/max in milliseconds:

| Atoms | Active image pairs | Brute force | Cells | Median allocated bytes, brute / cells |
|---|---:|---|---|---|
| 32 | 16 | 0.0039 / 0.0042 / 0.0051 | 0.0115 / 0.0125 / 0.0175 | 12,160 / 33,184 |
| 256 | 631 | 0.0782 / 0.0791 / 0.3675 | 0.1209 / 0.1256 / 0.1659 | 74,904 / 353,360 |

Raw printed report: `/tmp/pfw-choose-example.log`. These are diagnostic geometry
timings for the deterministic fixtures, not complete force-evaluation speedups
or evidence for a universal threshold. No default threshold was changed.

## NeighbourLists-style query interoperability (2026-10-21)

Inspected upstream README/current API and locally installed NeighbourLists 0.6.2
cell_list.jl, gpu_kernels.jl and AtomsBase extension, plus 0.5.10 legacy query
methods. The raw API stores lattice vectors in ROWS and uses
`R = X[j] - X[i] + transpose(cell)*S`. Current `neighbours` returns THREE arrays
(indices, displacements, shifts); 0.5's own `neighbours` historically returns two,
with its three-array query named `neigss`. This adapter targets the current
three-array contract, including on Julia 1.10. No upstream implementation was
copied. Source: https://github.com/JuliaMolSim/NeighbourLists.jl .

Added `neighbour_list(X, cutoff, cell, pbc)` and `neighbours(list,i)` as generics
owned by NeighbourList and reexported at the package root. NeighbourInterface
implements this boundary after Execution is loaded, reusing prepare!, canonical
cell topology and SIMD geometry rather than adding an all-pairs implementation.
Raw vectors of three-vectors and 3×N matrices are accepted; a small internal
AtomsBase adapter bridges numeric input to the existing system boundary. The
raw row-vector cell is transposed once. An AtomsBase system + cutoff overload
is also supplied. Numeric lengths mean Å, and unitful lengths convert explicitly;
query displacements are always numeric Å. Float64/width 4/Int32 indices and
shifts are defaults, with precision/width/bruteforce_threshold/int_type options.
Signed integer narrowing is checked, including inverse-shift overflow.

Queries use original atom indices and materialize a directed geometry snapshot.
Each owned half-list lane creates both directed endpoints, including ±self
images. Zero-shift self pairs remain excluded. Topology CellImage entries now
retain exact integer lattice shifts alongside Cartesian shifts (24 additional
bytes per candidate cell-image entry on this 64-bit machine). Caller-coordinate
windings are computed using the identical internal basis and precision as
prepare!, so reported shifts refer to original, possibly unwrapped positions.
Minimum-image corrections are explicitly included for the orthorhombic fast
path. This avoids reconstructing integer images from approximate displacement
vectors or silently reporting shifts of wrapped coordinates instead of inputs.
Open-axis shifts remain zero, including zero/infinite vacuum vectors.

Snapshot storage consists of per-original-atom vectors of indices, displacement
static vectors and integer shifts. Geometry/inputs are copied at construction;
queries return independent arrays with bounds checks and no ordering guarantee.
It is separate from the evaluation Workspace, and callers rebuild after changes.
This materialized compatibility boundary uses O(directed pair images) geometry
storage; ordinary force evaluation continues immediate block accumulation without
retaining pair results or allocating per-thread forces. No lazy-query/GPU backend
or automatic mutation/rebuild interface is claimed.

NeighbourLists is an OPTIONAL weak dependency, compat 0.5.10/0.6, with an extension
implementing its neighbours/nsites/npairs methods ONLY for our snapshot type.
External upstream constructor methods and upstream types are unchanged; there is
no method piracy. Three-array results follow our current-API contract even when
the loaded upstream version is 0.5. The dependency is a direct test extra for
comparisons; package loading and our own query interface need no NeighbourLists
installation. Julia compatibility/package identity are preserved. README shows
row-cell conventions, units, snapshot lifetime, qualified usage and limitations.

Validation compares directed images to the independent exhaustive oracle over
orthorhombic, general skewed, tiny multi-image and partially/nonperiodic cells,
checking original-coordinate reconstruction, exact integer labels, opposite
edges, ownership, brute/cell and scalar/SIMD consistency. Additional checks
compare directly with upstream query results for a nonsymmetric skew cell,
verify extension dispatch, bounds/errors, empty/single atoms, units, Float32,
custom integer types/overflow, snapshot independence and self-image directions.

Reference-version nuance observed during validation: 0.6.2 matches original-input
integer image labels directly. Legacy 0.5.10 wraps positions internally and its
`neigss` labels refer to those canonical positions; comparing them unchanged to
original-position shifts fails while displacement vectors agree. The legacy
reference test therefore supplies canonical periodic positions and converts its
labels with the fixture's known integer windings. Independent oracle and
original-coordinate reconstruction assertions remain unchanged on both versions.
This documents the upstream semantic difference rather than changing our output
or dropping the unwrapped-position checks.

Final suites pass **14,970 assertions** per process on Julia **1.13.2 / 4 threads**
with NeighbourLists **0.6.2**, and Julia **1.10.12 / 4 threads** with NeighbourLists
**0.5.10** (9,960 added assertions). AMD Ryzen 9 9950X, 16 cores/32 logical CPUs;
AtomsBase 0.5.3, AtomsCalculators 0.2.3, SIMD 3.8.1, StaticArrays 1.9.23,
Unitful 1.29.0, ForwardDiff 1.4.6, ChainRulesCore 1.26.4, Zygote 0.7.13 and
EmpiricalPotentials 0.2.6. Julia 1.10 uses the isolated source/environment described
above; root Julia 1.13 Manifest is preserved. Logs:
`/tmp/pfw-interface-113-final.log`, `/tmp/pfw-interface-110-final.log`.
A standalone root-environment process also passed the documented example with
the NeighbourLists extension absent, verifying the optional base interface.
Image-winding arithmetic and reverse-shift negation use checked integer
operations in addition to checked narrowing; an Int8 -128/+128 regression checks
inverse-image overflow. `git diff --check` passes. No performance measurements
were added or prior speedups attributed to this materialized query adapter.

## Two functional source submodules (2026-10-21)

At the user's request, the current source layout supersedes the earlier helper
module boundaries recorded above. `src/PotentialFrameWork.jl` now includes and
reexports only `NeighbourList` and `PairPotentials`, with their entrypoint files
under `src/submodules/NeighbourList/` and `src/submodules/PairPotentials/`.
The singular NeighbourList name and package identity are preserved.

NeighbourList includes storage.jl (workspace, cells, topology and scheduling),
geometry.jl (SIMD block packing and scalar/SIMD select), interface.jl (directed
snapshot queries) and benchmarks.jl (strategy comparison). PairPotentials imports
these neighbour helpers and includes models.jl (callbacks/adapters), evaluation.jl
(fused accumulation, prepared execution and public evaluation), calculators.jl
(AtomsCalculators boundary/workspace selection) and parameters.jl (energy
parameters and derivative hooks). The dependency is one-way from PairPotentials
to NeighbourList. SIMD masking is defined once in NeighbourList and shared.

Execution, CalculatorInterface, ParameterEnergy, NeighbourInterface and
NeighbourBenchmarks no longer exist as package submodules; they are not retained
as aliases. Their contents live directly in the two functional modules. Public
package-level evaluation, neighbour query, benchmarking and parameter APIs retain
the same names/signatures. Internal geometry paths move to NeighbourList;
cell_task!/execute!/outputs and parameter_seed move to PairPotentials. Tests,
benchmark/code-inspection scripts and current benchmark documentation use these
new locations. Historical reports above retain their original module names.

Julia package-extension modules remain in ext/, as required by Julia's extension
mechanism, and now import snapshot/parameter helpers from their functional owners.
They do not introduce additional modules beneath the package's src hierarchy.
The root module is limited to includes, explicit imports and exports. README
lists the current file tree and responsibility split. No numerical algorithms,
list ownership, unit conventions, topology lifetime or accumulation strategy
were changed. No new dependency or API compatibility requirement was introduced.

A standalone load confirmed the removed module bindings are absent and that
neighbour_list/evaluate! belong to NeighbourList/PairPotentials respectively.
Full numerical, scheduling, interoperability and AD regressions are used for
validation rather than adding tests that mirror the directory layout.

Final Pkg.test results: **14,970 assertions passed** on Julia **1.13.2 / 4 threads**
and **1.10.12 / 4 threads**, on AMD Ryzen 9 9950X (16 cores/32 logical CPUs).
The suites include numerical/virial/force derivatives, SIMD widths/precisions,
phase scheduling/write sets, directed image interoperability, strategy benchmark
validation, AtomsCalculators and ForwardDiff/ChainRules/Zygote regressions.
Dependency versions are unchanged from the preceding interoperability validation,
including NeighbourLists 0.6.2 on Julia 1.13 and 0.5.10 on Julia 1.10. Logs:
`/tmp/pfw-layout-113-t4.log`, `/tmp/pfw-layout-110-t4.log`. The Julia 1.10 source
snapshot is `/tmp/pfw-layout-110`, with its compatible isolated Manifest; the
root Manifest was preserved. Final `git diff --check` passes and active source,
examples, tests, extensions and benchmark scripts contain no stale helper-module
references. Historical generated benchmark/code-inspection artifacts were
preserved; no new performance or code-generation result is claimed.
