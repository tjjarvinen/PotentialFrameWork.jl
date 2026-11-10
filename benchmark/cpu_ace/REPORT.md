# CPU ACE optimization stage 1 — 10 November 2026

Implemented unrestricted prepared-CPU products using a division-free prefix and
reverse sweep. `product_algorithm=:reference` preserves polynomial projection;
`:auto` and `:optimized` use the fast path only for `orders === nothing`.
Nonempty U/V partitions, repeated indices, zero factors and current weights are
supported. Restricted policies always use the general algorithm. Functional AD,
equivariant features, model imports/checkpoints and GPU mathematics retain their
existing paths. Production CPU force ownership and scratch scheduling are unchanged.

## Source, fixtures and validation

The tracked starting tree was clean at
`ea076d3c18d5ae66776f1fa785e729f2eadfba54`; only the supplied instructions were
untracked. Its immutable archive was captured before editing. See
[baseline reproduction](baselines/README.md), [source patch](results/2026-11-10/implementation.patch),
[initial harness snapshot](results/2026-11-10/harness-and-new-tests.tar.gz),
[environment snapshots](results/2026-11-10/environments), and
[original input hashes](results/2026-11-10/original-input-sha256.txt).
Historical benchmark files and original `data/ace_ref` inputs were preserved.
Historical-source timings below were collected in this experiment, using the
same dependencies as the modified checkout; they are not relabelled old timings.

Modern references use ACEpotentials **0.10.2**, whose `ace1_model` constructs
`Models.ACEModel` with `EquivariantTensors.SparseACEbasis`, not old ACE1.
Medium and large models have deterministic nonzero species-dependent synthetic
weights and no physical-accuracy claim. Every product has an active precontracted
weight for every central species in these fixtures. The supplied small fitted
TiAl model is retained without refitting.

| Modern model | Species | Radial / pair | Harmonics | Density | Product orders and counts | Coupling shape / nnz |
|---|---|---:|---:|---:|---|---|
| small TiAl | Ti, Al | 30 / 12 | 9 | 27 | 1:12, 2:63, 3:92 | 123×167 / 167 |
| medium Ti | Ti | 37 / 10 | 25 | 43 | 1:10, 2:81, 3:139 | 110×230 / 230 |
| large TiAl | Ti, Al | 108 / 24 | 49 | 165 | 1:24, 2:630, 3:4264, 4:9001 | 3727×13919 / 14855 |

All modern models use exact exported cubic spline/transform/envelope data at
5.5 Å cutoff. Supplementary published Ti and Ti–Al–V retain exact **analytic**
radials and use independent **ACEpotentials 0.6.6 / ACE1** oracles. Converted
product counts are 16,444 / 51,379, densities 329 / 441, coupling shapes
12539×16444 / 130755×51379, and nnz 30,017 / 271,722. Their PI coordinates
are not individual invariant fitting rows. Attribution and licensing are retained
in the [fixture cache](results/2026-11-10/models/README.md).

Both full `test/runtests.jl` suites pass **54,593 assertions**, on Julia **1.13.3**
and **1.10.12**, four threads. This includes **635 new assertions** and existing
force-loss/MD, parameter/position, equivariance/parity, Lux, lifetime, scheduling,
units, import/persistence and portable GPU-kernel regressions. Optional ACEfit
0.3.1 passes four assertions. A separate gate passes **44 public E/F/V cases**,
including Float32 large-model 11,664-atom predictions. New prepared force/strain
stencils retain the 2e-7 bound at h=1e-5. Product finite differences use tolerances
appropriate to O(h²)+O(eps/h), separately from tighter algebraic comparisons;
existing tolerances were unchanged. See [validation](results/2026-11-10/validation).
CPU-only loading passes with AMDGPU/SpheriCart extensions absent. Actual RX
hardware and optional UnicodePlots were not rerun; no new coverage is claimed.

## Measurements

AMD Ryzen 9 9950X, 16 physical cores / 32 logical CPUs; Julia 1.13.3,
LLVM 20.1.8, BenchmarkTools 1.8.0, SIMD 3.8.2, SpheriCart 0.2.4,
StaticArrays 1.9.24, KA 0.9.45 and AK 0.4.4. BLAS uses one thread. Expensive
cases ran sequentially. Initial setup-contended observations are explicitly
labelled exploratory; final observations were collected separately.

Prepared combined E/F/V excludes extraction, graph/executor preparation and
owning output conversion, and includes refreshing coefficients from current
parameters. Selected **serial Float64** medians:

| Model | Atoms | Reference ms | Optimized ms | Incremental gain | Historical-source gain |
|---|---:|---:|---:|---:|---:|
| small TiAl | 128 | 1.743 | 1.317 | 1.32× | 1.33× |
| small TiAl | 11,664 | 159.648 | 120.121 | 1.33× | 1.32× |
| medium Ti | 1,024 | 17.847 | 13.014 | 1.37× | 1.38× |
| large TiAl | 128 | 72.907 | 17.584 | 4.15× | 4.15× |
| large TiAl | 1,024 | 578.571 | 137.960 | 4.19× | 4.27× |
| large TiAl | 11,664 | 6671.578 | 1604.239 | 4.16× | unmeasured |
| published Ti | 128 | 74.664 | 24.805 | 3.01× | 3.00× |
| published Ti–Al–V | 128 | 217.118 | 60.118 | 3.61× | 3.61× |

This is a substantial incremental improvement, not a universal win over modern
ACEpotentials. At 11,664 atoms, complete public Float64 calls take **317.7 ms PFW
versus 232.7 ms modern** for small TiAl, and **1798.9 versus 1383.3 ms** for large
TiAl. With four threads pinned to physical CPUs 0–3, large TiAl at 1,024 atoms
instead takes **74.0 versus 90.6 ms**; at sixteen threads on CPUs 0–15 it takes
**56.1 versus 49.7 ms**. Both sides use matched models/precision and fresh public
neighbor preparation in these modern comparisons. See [controlled repeats](results/2026-11-10/pinned/REPORT.md).
Unpinned parallel historical/reference runs differ appreciably between processes;
placement/compiler/executor effects cannot all be attributed to products.

Every final prepared/public trial has ten warmed samples. The large-model
11,664-atom reference **three-frame replay has five**, explicitly flagged under
the 90-second cap. Replay includes three perturbed frames and fresh rebuilds on
every frame, with owning outputs; no CPU skin reuse was introduced. Legacy
external timings instead retain JuLIP's native fixed-geometry neighbor cache and
three separate E/F/V calls; they are a different boundary. Float32 external
execution is unsupported, explicitly labelled rather than silently promoted.

For large TiAl at 11,664 atoms, Float64 force disagreement with the independent
modern evaluator is about 4.3e-14 eV/Å. For small TiAl, Float32 application
sensitivity against Float64 reaches 1.7e-5 eV/Å and 0.055 eV in the extensive
virial at that size; the rounded-PFW reference/optimized comparison is separate.
Float32 does not deliver a general 2× prepared-evaluator gain. Width experiments
at 128 atoms favor width 1 here: small TiAl optimized Float64 takes 1.223 / 1.317 /
1.818 ms at widths 1 / 4 / 8. No global width default was changed.

Serial prepared trials allocate **zero host bytes** in these measurements.
Parallel scheduling still allocates tasks: large TiAl / 1,024 atoms allocates
456,864 bytes at four threads and 1,002,032 at sixteen. Scratch payloads are
888,256 / 3,553,024 bytes respectively, plus one shared 24,576-byte force array.
These are allocation volumes and array payloads, not allocator peaks.

## Diagnostics and next stage

Isolated warmed first-site product/reverse time falls from approximately
558 to 98 µs on the large modern fixture. Its radial, density and force
back-contraction stages measure roughly 11.3 / 5.9 / 5.2 µs. Published Ti–Al–V
still spends about 370 µs in products/reverse versus 12.4 / 41.8 / 29.0 µs in
those stages. Sampling profiles support investigating products next for these
models. Stage measurements perturb caches and must not be summed into a total.
Small-model inference instead has substantial radial and full preparation costs;
parallel numerical speedups expose those costs further.

Actual Float32/Float64 radial and product callee dumps show precision-preserving
arithmetic. The prefix reverse is division-free, scalar and indexed, with bounds
checks retained; steady serial allocation measurements are zero. Dispatcher dumps
are distinguished from actual callee dumps. Static instructions are not runtime
percentages. General feature workloads include scalar, both vector parities and
rank-two tensor channels, plus existing joint position/parameter pullbacks. At
32 atoms, serial feature evaluation takes 0.146 ms and pullback 2.65 ms; these
allocating functional paths are preservation baselines, not accelerated AD.

Prioritize **shared-product plans/DAGs** for product-heavy models, retaining
original parameter coordinates and validating derivatives separately if shared
with training. Inspect structural species plans, without pruning fitted zeros
or assuming inactive product weights (all are active in these fixtures). For
small-model public calls, profile species packing/snapshot preparation and
scheduling separately before broader executor work. Across-site SIMD and broad
contraction/reduction redesign remain deferred hypotheses.

PACE is available and its executable/help/version/hash were probed. No timing
comparison was launched: equivalent exported spline interpolation, endpoints,
pair tables and E/F/V conventions are not validated. This is a structured optional
unsupported result, not a passing numerical comparison.

Reproduce with [README commands](README.md) and
[exact run commands](results/2026-11-10/REPRODUCE.md). Raw observations and expanded
tables: [quick](results/2026-11-10/final-quick/REPORT.md),
[1,024/thread scaling](results/2026-11-10/full-1024/REPORT.md),
[large structures](results/2026-11-10/large-11664/REPORT.md),
[legacy](results/2026-11-10/legacy-128/REPORT.md),
[widths](results/2026-11-10/widths/REPORT.md).
