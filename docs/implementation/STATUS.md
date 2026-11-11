# Implementation status

## CPU ACE product DAG (2026-11-10)

Completed the supplied `CPU-ACE-Product-DAG-Codex-Instructions.md` scope from
`4985499b6eb12421f512c3cb9d27671e5f76ec9b`. Tracked starting source was clean;
the supplied untracked instructions were preserved with the immutable source
archive identity. Historical stage-one results, reference inputs, model/oracle
caches and pinned manifests are unchanged. No agents were delegated.

Prepared scalar CPU ACE now accepts `product_algorithm=:dag` alongside preserved
`:reference`, `:optimized` and `:auto` behavior. A deterministic canonical-prefix
plan shares scalar products and reverses weighted terminal seeds additively,
preserving every original coordinate, repeated factors, zero factors and weights.
One parameter-independent plan is shared read-only; numerical node values and
adjoints belong to exclusive scratch slots. Coefficients/masks refresh weights
without rebuilding the plan; structural snapshots reject stale product topology.
CPUExecutor and generic NeighbourInterface preparation, growth and legacy
constructors are covered. Restricted-order requests retain polynomial fallback;
`ACE.effective_product_algorithm` reports the prepared effective path.

Auto stays on prefix/reverse. Functional AD and equivariant features retain
polynomial projection and their existing derivative rules; this mutation path
adds no AD claim. Radial/SIMD defaults, neighborhoods, scheduler, one shared force
array, units, model/checkpoint coordinates and GPU execution are unchanged.

Validation: **56,658 assertions / 62 testsets pass per full process**, Julia
**1.13.3 and 1.10.12**, four threads (2,065 added assertions). Focused DAG checks
pass 2,700 assertions; optional ACEfit 0.3.1 passes four, and CPU-only loading
passes with vendor/angular/plot extensions absent. Force/strain differences keep
2e-7 tolerance at h=1e-5; product/weighted-adjoint tolerances are unchanged.
Repeated inputs/weights/algorithms, zeros/repeats/shared descendants, nonzero U/V,
pruning, custom radials, growth, constructor compatibility, center selection,
generic batching and existing AD/geometry/lifetime/equivariance regressions pass.
Direct full-model comparisons against immutable stage-one output pass **72 cases**
across all fixtures, both precisions, 128/1,024 and selected 11,664-atom cases.
Existing prefix outputs are identical. DAG maximum force/virial disagreement is
2.504e-14 / 1.819e-11 in Float64 and 1.336e-5 / 0.00390625 in Float32, in eV/Å
and eV respectively. Independent modern/legacy cached oracles pass where available;
legacy independent caches cover 128 atoms only. Explicit arbitrary-order
reassociation errors and Float32 application sensitivities are recorded.

On Ryzen 9 9950X, Julia 1.13.3, pinned cores and BLAS=1, **612 sequential
BenchmarkTools trials all reached ten warmed samples**. The grid covers all five
fixtures at 128/1,024 atoms in both precisions, small/large modern 11,664-atom
cases in both precisions, real-site product/forward/reverse/readout diagnostics,
prepared/public/fresh-calculator calls, immutable stage-one comparisons and
1/4/16-core scaling. Width 4, fusion off and SpheriCart are fixed.

Large TiAl Float64 at 1,024 atoms improves prefix → DAG **139.909 → 68.453 ms
prepared (2.04×)** and **155.992 → 83.257 ms public (1.87×)**; at 11,664,
**1,663.012 → 775.673 ms prepared (2.14×)** and **1,845.275 → 975.139 ms public
(1.89×)**. Small/medium prepared gains are about 5–7%; small TiAl Float64 public
at 11,664 is effectively unchanged. Small-model parallel public medians slightly
regress within observed variation. Large TiAl 1,024 prepared gains at 1/4/16 cores
are 2.04× / 1.54× / 1.83×; small parallel cases are noisy. Shared integration raises
fresh prefix-calculator costs versus immutable stage one; these are reported,
not hidden in same-checkout speedups. No universal Float32/parallel gain is claimed.

Serial warmed numerical execution allocates zero host bytes in every measured
main/large case. Large TiAl has 16,229 multiply nodes / 13,919 output coordinates;
plan-array payload is 771,656 bytes, Float64 DAG scratch 262,304 bytes per slot,
retained polynomial tables 400 bytes and total slot scratch 484,368 bytes at
128/1,024 atoms. Sixteen slots use 7,749,888 scratch bytes plus the single 24,576-byte
shared force array at 1,024 atoms; parallel scheduling still allocates. Plan build
is 1.149 ms / 7,333,584 allocated bytes. All model/type counts and constructor costs
are retained separately from host allocation volume and allocator peaks.
Concrete forward/reverse callees show no observed Float32 widening or direct
allocator calls, and six warmed callee probes allocate zero bytes. Index checks
remain. `perf` is absent; no cache/bandwidth saturation claim is made. Portable
KA CPU regression tests pass; real GPUs and optional UnicodePlots were not rerun.

Report and review recommendation: `benchmark/cpu_ace/DAG-REPORT.md` (linked from
benchmark README). Exact commands, manifests, source/model/geometry hashes, raw
TOML/JSON, immutable output comparisons, reassociation evidence, profiles and code
inspection are under `benchmark/cpu_ace/results/2026-11-10-dag/`; reproduce with
`commands.sh` then `complete-validation.sh`, and the two full-suite commands in
`REPRODUCE.md`. The measured next review is large-model DAG reverse/indexing;
small-model public preparation is a distinct issue. Scalar sharing does not
establish MACE tensor, training AD or GPU speedups. No subsequent optimization
stage was implemented.

## CPU ACE optimization stage 1 (2026-11-10)

Completed the stage-one implementation from the supplied instructions on
`ea076d3c18d5ae66776f1fa785e729f2eadfba54`. Before editing, archived the clean
tracked source; the supplied untracked instruction file and original reference
inputs remain untouched. Immutable historical source, compatible application
manifests, source hashes/patches, cached model/geometry identities and exact
commands are recorded under `benchmark/cpu_ace/`. Historical benchmark artifacts
are unchanged. No delegation was used; package identity and Julia 1.10.10 core
compatibility remain unchanged.

Prepared CPU scalar inference now accepts `product_algorithm=:reference`,
`:optimized`, or `:auto` in executor/calculator construction, execution and generic
potential preparation. The reference uses the original polynomial algorithm.
The other choices use division-free prefix/reverse evaluation only when
`policy.orders === nothing`, including nonempty U/V species partitions.
Repeated/zero factors and weighted additive U/V adjoints use existing bounded
scratch, without per-product allocation or polynomial-table clearing. Restricted
policies always retain general projection. Parameters/masks refresh every call.
Functional AD, equivariant features, original imports/checkpoints, GPU mathematics,
CPU coloring, exclusive scratch ownership and one shared sorted force array are
preserved. The explicitly typed original five-argument PreparedACE constructor
remains compatible.

CPU-only quick/full BenchmarkTools entry points now cover matched modern
ACEpotentials 0.10.2 models: supplied small TiAl, deterministic synthetic medium
Ti and large TiAl. Their radial/harmonic/density/product counts and coupling
structure are recorded, including the large model's 13,919 products and 14,855
coupling nonzeros. Its actual backend is Models.ACEModel / EquivariantTensors,
not old ACE1. Published Ti/Ti–Al–V remain separately labelled 0.6.6/ACE1 analytic
comparisons, with attribution and exact pair/reference contributions retained.
The legacy environment runs on Julia 1.10.12; its bundled MbedTLS_jll manifest
cannot be loaded directly on Julia 1.13. Original pinned environments are intact.

Final sequential measurements cover 128/1,024/11,664 atoms, both PFW precisions,
1/4/16-thread processes, additional controlled physical-core affinity repeats,
widths 1/4/8, prepared/public/three-frame fresh-rebuild boundaries, stage
microbenchmarks, sampling profiles and actual Float32/Float64 callee dumps.
All final prepared/public trials have ten warmed samples; the expensive large
11,664-atom reference replay has five, explicitly flagged. Medium/large Float32
at 11,664 atoms are numerically validated but not timed. Legacy timings are
128 atoms only. No CPU skin reuse, product DAG, broad batching/contraction/
reduction rewrite or AD optimization was introduced.

On AMD Ryzen 9 9950X / Julia 1.13.3 / one BLAS thread, serial Float64 large TiAl
prepared E/F/V improves **578.57 -> 137.96 ms at 1,024 atoms (4.19×)** and
**6.672 -> 1.604 s at 11,664 (4.16×)**. The historical-source 1,024-atom result is
589.49 ms, a 4.27× total gain. Small TiAl improves **159.65 -> 120.12 ms** at
11,664 atoms (1.33× incremental, 1.32× historical). Published Ti and Ti–Al–V
improve 3.01× / 3.61× at 128 atoms. Complete public calls do not universally beat
modern ACEpotentials; losses and affinity-dependent parallel variation are
reported. External Float32 is unsupported and labelled separately from PFW
rounded-reference agreement and Float64 application sensitivity.

Serial prepared execution allocates zero measured host bytes in these cases;
parallel task scheduling still allocates (large TiAl / 1,024 atoms: 456,864 bytes
at four threads, 1,002,032 at sixteen). Scratch payload is separate from allocation
volume and allocator peaks. Product work remains dominant in large fixtures;
shared-product plans are the measured next priority. For small-model public
calls, preparation/species packing and scheduling merit separate profiling.
There is no universal Float32 or width speed claim. General scalar/vector/tensor
feature and position/parameter pullback costs are preserved as separate baselines.

Validation: **54,593 assertions pass per full test/runtests.jl process**, Julia
**1.13.3** and **1.10.12**, four threads (635 added assertions). Required existing
energy/force/virial derivatives, force pullbacks, fixed-cell rollout, equivariance/
parity, parameter authority, lifetime, units and portable GPU regressions pass.
The new prepared force/strain stencil retains 2e-7 tolerance at h=1e-5; product
finite-difference tolerances are separately justified, with existing tolerances
unchanged. Optional ACEfit 0.3.1 passes four assertions. Separate **44 public
E/F/V cases** pass against independent oracles and prepared outputs. CPU-only
loading passes with vendor/angular extensions absent. Real RX hardware and
optional UnicodePlots were not rerun; no new coverage is claimed.

PACE at `/opt/lammps-static/bin/lmp` is available and probed, with executable
hash/version/help retained. No PACE timing comparison was launched: equivalent
exported spline interpolation, endpoint/pair-table behavior and numerical E/F/V
conventions remain unvalidated. This is an explicit optional unsupported result.

Report, measured limitations and next-stage recommendation:
`benchmark/cpu_ace/REPORT.md`. Reproduction:
`benchmark/cpu_ace/README.md` and
`benchmark/cpu_ace/results/2026-11-10/REPRODUCE.md`. Raw TOML/JSON, source snapshots,
environments, validation logs, profiles, compressed generated code and host
identity are under `benchmark/cpu_ace/results/2026-11-10/`.


## GPU ACE milestone 1 (2026-11-09)

Implemented scalar GPU ACE E/F/V, coherent Float32/Float64 model/parameter
adaptation, optional AMD device integration, resident preparation and owning or
borrowed outputs, site batching, analytic callable radials, precision reports,
public fitted-model conversion/oracles and a sequential BenchmarkTools workflow.
The root package identity, Julia 1.10.10 compatibility declaration, existing CPU
AD/equivariant paths, local pure protocol and pinned reference data are preserved.
The supplied milestone instruction file remains untouched.

The RX 6800 XT is available: the restricted sandbox hid device nodes. Actual
hardware checks use the host-accessible device, not a CPU kernel substitute.
The tested stack is Julia 1.13.3, AMDGPU 2.10.0, KernelAbstractions 0.10.7,
AcceleratedKernels 0.5.12 and HIP 7.1.52802 on gfx1030 (72 CUs). CPU checks also
use Julia 1.10.12; the independent SpheriCart 0.2.4 checks use a separate KA 0.9
environment. No mandatory vendor dependency was added.

Validation: full suites pass **53,958 assertions on each Julia version**; actual
RX hardware passes **758 ACE assertions** and **1,539 existing pair assertions**.
The public Ti and April08 Ti–Al–V models pass **34 independent checks per backend**
(CPU and AMD), including both precisions and small/all-species and 128-atom bulk
fixtures. Angular values/Cartesian derivatives include independent SpheriCart
comparisons at lmax 0/1/2/4/6. Force and strain finite differences at h=1e-5
have worst errors 1.919e-7 eV/Å, 1.587e-7 eV (open strain), and 5.157e-8 eV
(periodic diagonal strain), against the unchanged 2e-7 absolute tolerance.
Net force norm is 7.217e-16 eV/Å. Full-order projection gives zero force difference.

Public legacy JSONs contain precontracted complex PI coordinates, not the paper's
independent invariant fitting rows. The isolated ACEpotentials 0.6.6 exporter
preserves fitted coefficients/order, pair/reference terms and exact analytic
radial recurrences, and expands complex products algebraically into real
harmonics. No refitting or tabulation approximation is used. Individual-row
invariant inspection/selection and invariant fitting reject this representation;
scalar inference, central selection and neighbor-species order projection work.
Source DOI, license, Crown/AWE notice, fit selection and checksums are retained in
`benchmark/gpu_ace/THIRD_PARTY.md` and acquisition metadata. TiAl's pinned original
ACEpotentials 0.10.2 model remains a separate oracle.

Both quick and bounded full presets ran on CPU and AMD. Quick includes all three
models at 128 atoms and TiAl at 11,664 atoms; full uses `--sizes 4` for both
precisions, denser cells, separate thread-count processes, staged/baseline/compact
and deterministic five-frame sequences. The larger default public-model and
27,648-atom cases remain unmeasured. In the quick Float64 run, staged prepared
TiAl at 11,664 atoms measured 29.44 ms AMD / 50.17 ms CPU (four threads). At
128 atoms, staged AMD is slower than CPU for all models. Compact reduces the
128-atom AMD prepared medians to 1.43 ms TiAl, 45.14 ms Ti and 251.56 ms alloy;
these short-budget timings include slow cases with one sample. Original-package
oracle boundaries and raw samples are retained; no universal speedup is claimed.

Automatic compact selection is restricted to the tested RX 6800 XT / Float64 /
128–1,024-atom small-spline range (lmax ≤ 2, ≤27 densities, ≤167 products).
TiAl at 1,024 atoms measured 11.07 ms staged / 3.845 ms compact; 11,664 atoms
show no robust compact advantage. Other configurations retain staged defaults
and manual overrides. The final automatic-path checks verify all three sizes.
CPU SIMD source and actual LLVM/native kernels are documented in
`docs/gpu_ace_cpu_simd_audit.md`: packed spline Horner and compiler-vectorized
SpheriCart coexist with mainly scalar density/product/reverse/scatter loops.
Batched SpheriCart remains the measured CPU baseline. Initial wrapper-only dumps
are preserved and explicitly distinguished from corrected inspection artifacts.

Reproduction and timing boundaries: `benchmark/gpu_ace/README.md`; public APIs and
limitations: `docs/gpu_ace.md`; raw evidence is under
`benchmark/gpu_ace/results/2026-11-09-julia1.13.3-5eef1c1-dirty/`.
GPU execution remains orthorhombic single-image scalar inference. GPU AD/training,
equivariant outputs, skew/multiple images and other vendor validation are deferred.
Device allocator high-water marks/register spills are not measured; reported
array payloads and host allocations are separate. CPU optional integrations not
exercised by the recorded full suites are not newly claimed as validated.


## CPU equivariant ACE feature layer (2026-11-04)

Implemented the feature phase from `ACE-Equivariant-Codex-Instructions.md` on
`13f5cf5e249692bbe0cdb0890d1476b09cdaf339`. The supplied untracked instruction
file, local reference data/environments and unrelated files are preserved.
`ACEFeatureLayer`, `feature_parameters`, `feature_basis_summary` and
`ace_features(X, layer, ps, environment, H)` live in ACE. The root module,
scalar ACE evaluator/prepared executor, calculators, import/checkpoint formats,
ACEfit path, CPU scheduling and GPU code are unchanged. No dependency or Julia
compatibility requirement was added.

The generated constructor accepts species, an explicit radial basis, correlation,
input lmax, totaldegree truncation and ordered `(L, parity, channels)` blocks.
It generates complete species-resolved density channels, commuting products and
sparse real Clebsch–Gordan multiplets, removes zero/dependent coupling paths as
whole multiplets, and reports available basis multiplicities separately from
readout width/rank bounds. Intermediate/output L may exceed input lmax. A small
resolved-table constructor and an identical-coordinate scalar 0e adapter provide
independent reference checks. The scalar adapter excludes reference/pair terms
and rejects species-order projected bases. Mixing parameters have an explicit
central-species axis and share weights across magnetic components.

Owning tuple outputs have `(2L+1, channels, original atoms)` axes. The serial
feature path reuses full directed environments, radial/angular evaluation and the
scalar division-free product reverse. The ChainRules rule contracts arbitrary
feature cotangents into position, mixing, radial-scale and callable-radial
parameter gradients without full Jacobians. It snapshots model/inputs/parameters/
topology for pullback lifetime safety. Zero/thunked blocks and empty/isolated
systems are covered. Custom parameter VJPs use ForwardDiff on locally contracted
radial values. Analytic position derivatives need the callable d/dr contract;
trainable callable values must declare `ad=true` and support generic ForwardDiff
arithmetic. Lux setup/apply uses the existing `lux_layer` convention with owning
explicit parameters and unchanged state.

Focused suites pass **163 assertions** on both Julia **1.13.3** and **1.10.12**,
four threads. They cover a hand-resolved nonzero cross-product 1e fixture,
generated 0e/1o/1e blocks, symmetric-cube multiplicity/dependence removal,
correlation three and L=3, real-harmonic rotations, inversion, translation,
permutation/original order, invalid/unavailable irreps, wider readouts, periodic
multi-image topology and self-images, arbitrary/single-block cotangents,
position/parameter finite differences and ForwardDiff checks, saved pullbacks
after later forwards and mutation, Lux, a composed invariant head's forces, and
scalar many-body equality with reference/pair subtraction.

The full suites on Julia **1.13.3** and **1.10.12** each pass **49,222 assertions**
(49,059 existing plus 163 feature assertions), including scalar ACE, import, persistence, fitting-row,
AD, neighborhood and portable CPU GPU-kernel regressions. The optional existing
scalar ACEfit regression passes **4 assertions** with ACEfit **0.3.1**. The
separate UnicodePlots tests were not rerun (UnicodePlots is not installed in the
local depot); no plotting code changed. A standalone root-environment import
passes with the SpheriCart extension absent, confirming optional CPU-only loading.
Logs: `/tmp/pfw-equivariant-full113.log`, `/tmp/pfw-equivariant-full110.log`,
`/tmp/pfw-feature-tests110.log`,
`/tmp/pfw-feature-optional-tests.log`, `/tmp/pfw-feature-cpu-load.log`.

Host: AMD Ryzen 9 9950X (16 cores/32 logical CPUs, znver5). Both application
environments resolve SpheriCart 0.2.4, ChainRulesCore 1.26.5, ForwardDiff 1.4.6,
LuxCore 1.5.3, Zygote 0.7.15, AtomsBase 0.5.3, AtomsCalculators 0.2.3,
StaticArrays 1.9.24, SIMD 3.8.1, Unitful 1.29.0, KernelAbstractions 0.9.45,
AcceleratedKernels 0.5.10 and NeighbourLists 0.5.10. The root and pinned reference
Project/Manifest files were preserved; tests used isolated application environments
and writable temporary depots layered over the installed packages.

`benchmark/ace_features.jl` measures the warmed serial reference path, with
preparation excluded. For its three-atom/two-species correlation-two fixture
(14 densities, 107 products, 15/12/6 coupled basis channels), Julia 1.13.3 / four
threads reports **16,816 bytes** for forward evaluation, **78,936 bytes** for
rrule forward including snapshot creation, and **79,088 bytes** for pullback.
Output array payload is **360 bytes** and saved-pullback `summarysize` is
**17,536 bytes**. These are allocating correctness paths; no speedup or
allocation-free claim is made. The diagnostic timings in
`/tmp/pfw-equivariant-costs.log` were collected while a separate regression run
was active and are not used as performance evidence. No new production force
buffers or CPU reduction strategy were introduced.

Reproduction with this checkout developed in an environment containing the test
extras (the application Project/Manifest files are under `/tmp/pfw-equivariant`
and `/tmp/pfw-equivariant110`):

```sh
JULIA_DEPOT_PATH=/tmp/pfw-julia-depot:/var/home/teemu/.julia julia --project=/tmp/pfw-equivariant -t4 test/ace_features_runtests.jl
JULIA_DEPOT_PATH=/tmp/pfw-julia-depot:/var/home/teemu/.julia julia --project=/tmp/pfw-equivariant -t4 test/runtests.jl
JULIA_DEPOT_PATH=/tmp/pfw-julia-depot110:/var/home/teemu/.julia /var/home/teemu/.julia/juliaup/julia-1.10.12+0.x64.linux.gnu/bin/julia --project=/tmp/pfw-equivariant110 -t4 test/runtests.jl
JULIA_DEPOT_PATH=/tmp/pfw-julia-depot:/var/home/teemu/.julia julia --project=/tmp/pfw-equivariant -t4 benchmark/ace_features.jl
JULIA_DEPOT_PATH=/tmp/pfw-julia-depot:/var/home/teemu/.julia PFW_ROOT="$PWD" julia --project=/tmp/pfw-equivariant-optional -t4 /tmp/pfw-feature-optional-tests.jl
```

See `docs/ace_features.md` for construction/evaluation signatures, parameter axes,
normalization/phase conventions, truncation and derivative contracts. Tables and
numerical reference computation use at least Float64 and retain promoted AD
arithmetic. Construction uses numerical multiplet rank detection and has
combinatorial cost; no large-basis optimization is claimed. Caller-supplied
resolved tables must themselves be valid intertwiners. Cell/topology derivatives,
feature SIMD/GPU, full MACE, general equivariant training and mixed derivatives
for composed-model force losses remain outside this phase. No new feature
checkpoint schema or scalar ACEfit replacement was introduced.


## GPU phase 1: sorted cells and radial pairs (2026-10-27)

Delivered G0–G5 on base `a6b3dd56e2560fe1339f0f41464da052f5cb7ff7`.
The user's appended AGENTS.md, untracked GPU specification/plan/validation files,
reference environments/data and package identity are preserved. GPU execution is
explicit opt-in inference; existing CPU APIs, multi-image/skew geometry, colored
force ownership, calculators and functional AD remain supported.

`NeighbourList` owns GPU configuration/workspace, generic geometry extraction,
sparse integer-triple record keys, AK record sorting, device unpack/ranges/binary
lookup, scalar NI views/outputs, both KA interaction kernels, reductions and skin
connectivity. `PairPotentials` owns GPU model preparation, runtime species/parameter
binding, SmoothPair and shifted Lennard-Jones mathematics, custom radial
certification and public owning-result boundaries. NeighbourInterface remains a
pure protocol; AtomsCalculators and the root module are unchanged. Core KA/AK
have no mandatory GPU vendor package; AMDGPU is confined to the application
validation environment. Julia compatibility remains 1.10.10.

The default is Float64, central-thread ownership and fresh cells on every call.
Explicit Float32 is supported. Experimental candidate tiling assigns one bounded
central tile to one workgroup, stages central positions, loads contiguous candidate
runs and reduces center forces in local memory. There are no inter-workgroup force
atomics, full force buffers per workgroup or dense distance matrices. Both paths
use the same unweighted radial callback and full directed force/half E/W convention.
Parameters refresh on every supplied call, retaining original model species indices
and full named/isotopic ChemicalSpecies identities. Numerical results own device
storage; convenience results own unitful original-order host arrays. Prepared
inputs borrow reusable upload storage and reject after later preparation on the
same workspace. Numerical execution reads one reduced scalar to reject singular
origins/nonfinite callback output, including force-only requests; it never reads
forces merely to satisfy the host convenience API. This synchronization is timed.

Optional skin was implemented after fresh traversal and tiling. Device count/scan/
fill builds exact conservative atom CSR candidates at cutoff+skin; reuse updates
positions in a fixed sorted permutation and filters at the physical cutoff.
Device maximum displacement plus a conservative threshold margin triggers rebuild.
Raw-coordinate differences conservatively handle wrapped jumps without inferring
winding. Metadata/species/cutoff changes invalidate connectivity; capacity grows
without truncation. Skin=0 remains default and no automatic tuning is provided.

Final full suites pass **49,059 assertions per process** on Julia **1.13.3** and
**1.10.12**, four threads, including all prior 47,518 CPU/interface/AD assertions
and 1,541 new kernel/interface/error/lifetime checks. Logs:
`/tmp/pfw-gpu-regression-113.log`, `/tmp/pfw-gpu-regression-110.log`.
These resolve KA 0.9.45/AK 0.4.4 with SpheriCart 0.2.4. Compatibility adapters
also support the KA 0.10.3/AK 0.5.10 device environment. Separate POCL CPU kernel
suites pass on both Julia versions; its original struct-reinterpretation compiler
failure was fixed by explicit species-field decoding. Standalone CPU-only import
passed without AMDGPU/CUDA loaded (`/tmp/pfw-gpu-cpu-load.log`).

The RX 6800 XT real-device suite passes **1,539 assertions** on Julia 1.13.3,
AMDGPU 2.10.0, KA 0.10.3, AK 0.5.10, GPUCompiler 2.18.0, GPUArrays 12.1.1 and
LLVM.jl 10.4.2. The device reports gfx1030, 15.984 GiB and wavefront 32. Host:
Ryzen 9 9950X, Linux 7.2.11-200.fc44.x86_64, HIP 7.1.52802. Relevant dependencies:
AtomsBase 0.5.3, AtomsCalculators 0.2.3, AtomsBuilder 0.2.3, AtomsSystems 0.1.2,
SIMD 3.8.1, StaticArrays 1.9.24 and Unitful 1.29.0. Diagnostics/test summaries
are persisted in `GPU-PAIR-HARDWARE.txt`; exact benchmark dependencies are in
`benchmark/gpu_pair_measurements.toml`. Tests include all property subsets,
precision overrides, independent converted-geometry sums, force/strain differences,
strict species accessors, parameter updates, origins, sparse/open/mixed/periodic
cells, aliases, tails/dense populations, alternate launches, skin trajectories and
invalidation, owning outputs and unsupported requests. GPU AD explicitly rejects.

Final synchronized warm Float32 argon full E/F/W medians (seven repetitions,
32 Julia threads for CPU comparisons): 10,976 atoms, central numerical **0.479 ms**,
public GPU **0.728 ms**, public CPU **4.741 ms**; candidate numerical **1.395 ms**.
97,556 atoms: central numerical **1.477 ms**, public GPU **2.906 ms**, public CPU
**28.655 ms**; candidate numerical **11.693 ms**. Numerical timings include fresh
cells/sort/search/reduction, scalar validation and owning device result copies,
excluding upload/host force copies. Public timings include extraction/upload and
unitful host outputs. All 48 cases, first-call/capacity costs, variation, stage
profiles and memory/host allocations are recorded. Candidate tiling is slower in
these fresh-cell samples. In the forty-step nonuniform 10,976-atom trajectory,
skin 0.5 rebuilt once: central amortized **0.256 ms** versus fresh **0.455 ms**;
candidate **0.242 ms** versus fresh **1.371 ms**. Rebuild/reuse paths were warmed
separately; every evolving step matches fresh GPU execution outside timing.
These measurements describe the tested cases, not general performance promises.

Reproduce from the repository (GPU application setup and callback examples are in
`docs/gpu_pairs.md`):

```sh
julia --project -t4 -e 'using Pkg; Pkg.test(;julia_args=["--threads=4"])'
julia --project=/tmp/pfw-gpu test/gpu_pairs.jl --amdgpu
julia --project=/tmp/pfw-gpu examples/gpu_pairs.jl
julia --project=/tmp/pfw-gpu benchmark/gpu_pairs.jl
julia --project=/tmp/pfw-gpu benchmark/gpu_pair_force_errors.jl
```

Julia 1.10 used `/var/home/teemu/.julia/juliaup/julia-1.10.12+0.x64.linux.gnu/bin/julia`
with isolated `/tmp/pfw-gpu110` and `Pkg.test("PotentialFrameWork";
julia_args=["--threads=4"])`; the separate pinned ACE reference environment was
not changed. No real-device gate remains open for the tested RX 6800 XT environment.
Earlier AMDGPU 2.8 attempts failed sort-kernel compilation on Julia 1.12/1.13;
a Julia 1.10 LLVM attempt crashed. Those GPU combinations are not claimed as
supported. GPU skew/multi-image, GPU AD, GPU ACE/MACE and other vendors remain
unsupported/untested as documented. The G1 radial correctness gate uses the
production device-built cells; no transferred-CPU-connectivity prototype is kept.
See `GPU-PAIR-RESULTS.md` and `docs/gpu_pairs.md` for tolerances, storage/lifetimes,
limitations and complete reproduction details. Historical measurements below are
unchanged.


## Distance-first pair callback refinement (2026-10-26)

At the user's request, pair callbacks now primarily consume `NI.distance(block)`
and write `NI.set_radial_derivative!(out,dEdr)`. `NeighbourList` forms
`g=(dEdr/r)*d` and performs force/virial reduction. `NI.displacement(block)` and
Cartesian `set_gradient!` remain optional; no mandatory trait or potential
supertype was added. The native accessor computes scalar/SIMD distances from safe
current block geometry. Cartesian legacy paths incur no mandatory square root
or extra stored distance array. Output setters select the radial or Cartesian
representation; inactive values remain unconsumed. The established regular-origin
contract remains: dE/dr=0 at r=0, or explicit rejection by the potential.

`RadialPotential` and `SmoothPair` prepare native radial adapters. SmoothPair's
explicit triangular parameter indices retain original model species coordinates,
even when provider packing uses a different species order. Updated explicit
parameters remain authoritative. Historical callable adapters, anisotropic pair
callbacks and the functional derivative implementation are retained. The portable
pair example and independent reference backend now use radial outputs, and a
separate distance-only test view deliberately has no displacement accessor.

New tests compare native/legacy radial forces and virials over Float32/64,
SIMD widths 1/4/8 and serial/threaded execution, check the regular origin and
explicit parameter updates, and retain the existing property-subset, anisotropic,
periodic-image and derivative regressions. Final validation and measurements
are recorded below; earlier interface timings remain historical evidence.

Final full suites pass **47,518 assertions per process** on Julia **1.13.3** and
**1.10.12**, each with four threads (`/tmp/pfw-radial-final-113.log`,
`/tmp/pfw-radial-final-110.log`). This adds 61 focused radial checks to the prior
47,457-check suite. The root-environment distance-first example also passes
(`/tmp/pfw-radial-example.log`). Existing force-loss, parameter/position AD,
anisotropic virial, units, scheduling and reference regressions remain passing.
Reproduce with the same isolated application environments/commands in the section
below, or run `test/neighbour_evaluation.jl` through the existing focused runner.
No new dependencies or reference files were modified.

Final warmed inference (`benchmark/inspect_neighbour_evaluation.jl`) reports
0 bytes for the native radial pair cell task, legacy Cartesian pair cell task and
ACE site callback; the radial pair `code_warntype` has `Body::Nothing` and no
`::Any` locals (`/tmp/pfw-radial-inference.log`,
`/tmp/pfw-radial-pair-warntype.txt`). Sequential three-atom prepared measurements
use the same hardware/dependencies/options as the interface measurements below
and are saved separately in `benchmark/radial_pair_measurements.toml`: legacy
Cartesian SmoothPair median **0.841 μs**, native radial SmoothPair **0.811 μs**,
and the uniform analytic RadialPotential adapter **0.686 μs**. Each allocates
**0 bytes per warmed call**. The adapter uses a simpler species-independent
function; these tiny measurements establish neither a general speedup nor a
throughput guarantee. No extra force array or stored distance array was added;
SmoothPair retains only a species-sized structural parameter-index map. Historical
interface measurements remain unchanged. Final `git diff --check` passes.


## CPU neighborhood evaluation interface (2026-10-26)

Implemented on `neighbor_list_interface`, starting from
`576862ac018957aaa4bde9640b098dca7d6fc0c9`. Preserved the user's modified AGENTS.md,
untracked neighborhood specification/plan/validation documents, package identity,
compatibility, optional dependency boundaries and all reference data. This section
supersedes earlier module-ownership descriptions only for the completed interface.

M0–M6 are complete. `NeighbourInterface` contains only declarations/contracts and
imports AtomsCalculators; it has no backend/model imports or allocation/traversal
algorithms. `NeighbourList` owns full connectivity/species packing, safe native
pair/site/batch views and supplied outputs, direct execution, loop batching,
colored phases and reduction. ACE retains mathematical channel/cutoff selection,
local scratch, radial/angular/projection kernels and authoritative parameter/cache
updates. `ACE.FullEnvironment` and `ACE.edge_displacement` remain compatibility
bindings. PairPotentials retains models/binding, legacy request/callback adapters,
public extraction/output and calculator boundaries; fused pair block traversal is
preserved. Backend loops contain no ACE type checks.

Final signatures, ownership and implementer contracts are in
`docs/neighbour_evaluation.md`; `examples/neighbour_evaluation.jl` implements pair
and angular-site callbacks without native view-field access or mandatory supertypes.
Preparation, scratch and optional derived-coefficient refresh hooks use explicit
parameters and declared units. All property subsets use canonical typed AC tuples
inside callbacks; legacy requests remain supported. Existing calculator result-key
order, state, additive destinations and public units are preserved. Unsupported
native unit conventions fail explicitly; convertible system input units still work.

`batch_size` is a positive integer construction option on ACE CPUExecutor,
ACECalculator and generic site preparation. Size 1 dispatches directly; larger
sizes dispatch to a specialized batch method or the CPU loop fallback. Serial
selected-center order, full neighborhoods, bounded exclusive scratch slots, distinct
pending outputs, tail handling, stale version/capacity/layout rejection and owning
original-order totals are tested. Parallel work batches within complete-endpoint
colored phases, retaining one shared sorted force array and barriers. No global
edge-gradient array or worker-sized full force arrays were added.

`test/neighbour_evaluation.jl` adds 2,700 checks: independent differently stored
O(N²) open-cell pair/angular-site evaluation, angular finite-difference force and
strain checks, explicit unsupported reference geometry, every requested subset,
duplicates/order, empty requests/environments/systems, precision, SIMD padding,
poisoned inactive output storage, batch 1/2/nondivisors/tails, direct vs specialized
batch dispatch, serial/threaded agreement, selected centers/order, parameter updates,
model transformations and layout lifetime. Existing image, scheduling, unit,
calculator, fitting/design, persistence, custom/source radial and AD tests remain.
The additional standalone `test/neighbour_saved_reference.jl` passes 597 checks
against all 33 saved full-precision ACEpotentials predictions with width 1/4 and
batch 1/2/5, including source/geometry hashes. It does not regenerate/refit data.

Validation evidence:

- Initial unchanged-source Pkg.test: 44,757 assertions on Julia 1.13.3, child
  process 32 threads (`/tmp/pfw-neighbour-baseline.log`). The initial attempt hit
  the read-only default depot; a writable `/tmp/pfw-depot` overlay fixed this.
  Registry DNS was unavailable, but installed packages plus offline resolution
  supplied all required core/test dependencies without modifying root manifests.
- Final Julia 1.13.3 / 4 threads: 47,457 assertions
  (`/tmp/pfw-neighbour-final2-113.log`). Final focused 1-thread suite: 2,700
  (`/tmp/pfw-neighbour-final2-focused-t1.log`).
- Final Julia 1.10.12 / 4 threads: 47,457 assertions
  (`/tmp/pfw-neighbour-final-110.log`), including optional NeighbourLists 0.5.10.
- Saved 33-frame check: 597 assertions (`/tmp/pfw-neighbour-saved-reference.log`).
- Dependency-light example passes in the root environment without loading any
  optional neighborhood/ACE/AD packages (`/tmp/pfw-neighbour-core-example.log`).
- `code_warntype` on the imported prepared site callback has `Body::Nothing` and
  no `::Any` locals; function-barrier warmed callback allocation is 0 bytes
  (`/tmp/pfw-neighbour-inference-final.log`, `/tmp/pfw-neighbour-site-warntype.txt`).
  Final pair cell-task inference also has `Body::Nothing`, no `::Any` locals
  and 0 warmed bytes (`/tmp/pfw-neighbour-final-inference.log`,
  `/tmp/pfw-neighbour-pair-warntype.txt`).
  The global-splat probe's 448 bytes were measurement-wrapper overhead, excluded
  from the function-barrier result. Legacy pair adapter inlining removed an
  initially observed 576-byte prepared-call allocation regression.

Development checks caught and fixed calculator output-key reordering, batch-method
fallback ambiguity, single-property result-shape handling and a dynamically built
prepared-evaluator type. Earlier test processes that loaded intermediate source
failed the later-added single-property tests; final clean processes use the fixed
source. Numerical tolerances and historical results were not weakened.

Hardware: AMD Ryzen 9 9950X, 16 cores/32 logical CPUs. Application environments use
AtomsBase 0.5.3, AtomsCalculators 0.2.3, SIMD 3.8.1, StaticArrays 1.9.24,
Unitful 1.29.0, SpheriCart 0.2.4, ForwardDiff 1.4.6, ChainRulesCore 1.26.4,
Zygote 0.7.14, LuxCore 1.5.3 and EmpiricalPotentials 0.2.6. NeighbourLists is
0.6.2 on Julia 1.13 and 0.5.10 on Julia 1.10. Julia 1.10.12 and 1.13.3 are the
available supported runtimes; exact 1.10.10/1.13.2 binaries are not installed.
The pinned Julia 1.13.2 reference environment is unchanged and was not regenerated.

Reproduction (isolated application environments include the Project test extras;
root and pinned manifests are preserved):

```sh
export JULIA_DEPOT_PATH=/tmp/pfw-depot:/var/home/teemu/.julia
julia --project=/tmp/pfw-neighbour-env --threads=4 test/runtests.jl
/var/home/teemu/.julia/juliaup/julia-1.10.12+0.x64.linux.gnu/bin/julia --project=/tmp/pfw-neighbour-env110 --threads=4 test/runtests.jl
julia --project=/tmp/pfw-neighbour-env --threads=4 test/neighbour_saved_reference.jl
julia --project --threads=1 examples/neighbour_evaluation.jl
julia --project=/tmp/pfw-neighbour-env benchmark/inspect_neighbour_evaluation.jl
julia --project=/tmp/pfw-neighbour-env --threads=4 benchmark/neighbour_evaluation.jl "$PWD" /tmp/pfw-neighbour-after.toml
```

For baseline timing, the unchanged source/project/manifest/ext snapshot is
`/tmp/pfw-neighbour-before`, with application environment
`/tmp/pfw-neighbour-before-env`. Activate that environment and pass the snapshot
as the benchmark's package-root argument. Run sequentially without concurrent tests.

Final sequential prepared measurements (Julia 1.13.3, four Julia threads,
BLAS one thread; serial execution; width 4, four exclusive slots; same three-atom
periodic fixture; 20 warmups then nine samples of 100 calls, excluding preparation
and owning outputs) are saved in `benchmark/neighbour_evaluation_measurements.toml`:

| Prepared evaluation | Before median | After median | Before/after bytes per call |
| --- | ---: | ---: | ---: |
| Pair E/F/W | 0.805 μs | 0.864 μs | 0 / 0 |
| ACE E/F/W, batch 1 | 35.416 μs | 35.780 μs | 96 / 0 |
| ACE E/F/W, batch 2 | — | 35.255 μs | — / 0 |
| ACE E/F/W, batch 5 | — | 35.105 μs | — / 0 |

The observed default-path differences are about +7.3% pair (59 ns) and +1.0%
ACE on this tiny fixture. Earlier repeated development samples varied enough
that these are not general slowdown/speedup estimates. No batch throughput gain
is claimed. Existing ACE scratch is 108,360 bytes across four slots before and
after. New batch-1 gradient storage is 2,304 bytes (24 padded edges/site capacity,
Float64, four slots); batch 2/5 require 4,608/11,520 bytes. Slot container summary
including the existing aliased scratch is 111,368 bytes at batch 1; do not add that
summary to scratch bytes as if they were disjoint allocations. The production
force storage remains one `3*N*sizeof(T)` array. Final `git diff --check` passes.
No required interface gate remains unverified in the available runtimes. Exact
pinned-binary re-execution and live ACEpotentials re-export were not performed;
the supplied saved oracle was checked without modification.

Remaining limits are intentional: mutating inference executors have no new AD
rules; the functional serial/ForwardDiff/ChainRules/Zygote/Lux path preserves the
existing E→X/ps, F→X/ps (including force-loss products), W→ps and documented W→X
capabilities/errors, saved inputs and repeat/delayed pullbacks. Legacy pair
adapters may compute/discard extras, and ACE prepared energy-only uses existing
joint derivative kernels. Scalar/custom radial fallback allocation limits remain.
Generic ACE preparation currently allocates a masked readout during its once-per-
execution refresh (the production CPUExecutor reuses its readout buffer), not per
site. Callback separation adds a bounded local gradient write/read and backend
geometry reconstruction. The batch loop is a correctness fallback with no optimized
contractions. No new site-energy API, skin reuse, gather/atomic reduction, GPU/MACE,
AD optimization, external protocol changes or upstream proposal was introduced.
The independent reference backend is deliberately limited to open geometry;
existing optional NeighbourLists query interoperability remains separate.


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
