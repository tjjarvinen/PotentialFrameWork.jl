# CPU SIMD audit — GPU ACE milestone 1

Measured 2026-11-09 on AMD Ryzen 9 9950X, `znver5`, Julia 1.13.3,
SIMD/SpheriCart/BenchmarkTools versions in the raw reports. Repository base
5eef1c1 with the milestone working tree; exact source hashes and dependency
snapshots accompany the results. These measurements preserve the historical
`benchmark/ace` artifacts.

All cases use the same deterministic 128-atom bcc workloads as the GPU suite,
serial prepared E/F/V, one Julia thread and one BLAS thread. Model/geometry
preparation and output copies are excluded. Each mode is warmed, `evals=1`,
up to five samples with a 0.12-second budget. Slow cases may have fewer samples;
raw samples are authoritative. The TiAl grid has 42 cases: widths 1/4/8,
scalar/SIMD, both precisions, SpheriCart/generated/block angular paths and valid
fusion alternatives. Each public analytic model has eight bounded width-4
cases covering scalar/SIMD, SpheriCart/block and both precisions.

| Model | Precision | Scalar radial + SpheriCart (ms) | SIMD radial + SpheriCart (ms) | SIMD radial + block (ms) | Host bytes, same column order |
| --- | --- | ---: | ---: | ---: | ---: |
| tial | Float64 | 1.846 | 1.758 | 1.889 | 3,088,384 / 0 / 288,768 |
| tial | Float32 | 1.764 | 1.663 | 1.784 | 1,916,928 / 0 / 151,552 |
| ti | Float64 | 74.892 | 74.798 | 78.625 | 7,364,608 / 0 / 1,080,320 |
| ti | Float32 | 73.404 | 73.706 | 74.902 | 4,751,360 / 0 / 547,840 |
| tialv | Float64 | 220.532 | 216.621 | 217.471 | 15,441,920 / 0 / 601,088 |
| tialv | Float32 | 215.248 | 213.594 | 218.393 | 8,314,880 / 0 / 317,440 |

The SIMD radial path avoids the allocating scalar `radial_ed` slices. For the
public callables, “SIMD radial” still evaluates each custom radial lane-wise;
it primarily supplies preallocated storage. The small time differences between
those modes do not establish vectorized custom mathematics. Float32 improves
some TiAl timings but provides little CPU benefit for these expanded legacy
products. There is no universal precision or speed claim.

Width-4 batched SpheriCart is the best measured TiAl mode in both precisions.
Generated angular and explicit blocks do not beat it here; width 8 and fusion
also fail to improve the best full E/F/V timing. The explicit block implementation
allocates recurrence arrays per call. Its vector arithmetic alone does not justify
promoting it. See every width/fusion case in `audit-corrected/tial/timings.toml`.

TiAl has 3,328 physical edges and 3,584 packed lanes (7.14% padding) at width 4;
the public models have 7,424 edges and 7,680 lanes (3.33% padding). TiAl's
72 colored phases contain one or two sites; both public workloads have 128
single-site phases. Thus extra executor threads expose almost no concurrency
for these small dense force-write footprints. The separate four-/one-thread
process results in the full preset agree with that scheduling constraint.

## Source and generated-code evidence

Optimized LLVM and native code are captured for concrete warmed signatures,
both precisions, for all three models. The initial audit's `spline` dump captured
a Julia default-argument wrapper. It is preserved in `full/audit-*`; **do not
interpret that wrapper as evidence of absent SIMD**. The corrected audit passes
explicit radial parameters and inspects site mathematics, force/virial reduction
and SpheriCart's actual normalization/harmonic/rescaling kernels directly.
The artifact stem `spline` names the generic radial inspection: for public models
its signature is a `CallableRadial`, not a spline.

The following counts are static occurrences of packed/scalar arithmetic mnemonics
in TiAl native text (`vadd`, `vmul`, `vdiv`, `vsqrt`, `vfmadd` families). They
are a compact index into the evidence, **not runtime instruction counts, SIMD
percentages or cycle attribution**. Calls and unlisted arithmetic are excluded.

| Inspected method | Precision | Packed occurrences | Scalar occurrences |
| --- | --- | ---: | ---: |
| spline | Float64 | 8 | 54 |
| angular | Float64 | 21 | 0 |
| density | Float64 | 0 | 12 |
| products_reverse | Float64 | 0 | 31 |
| site_math | Float64 | 19 | 66 |
| force_virial_reduction | Float64 | 0 | 13 |
| sphericart_normalize | Float64 | 22 | 13 |
| sphericart_kernel | Float64 | 64 | 75 |
| sphericart_rescale | Float64 | 12 | 5 |
| spline | Float32 | 8 | 54 |
| angular | Float32 | 21 | 0 |
| density | Float32 | 0 | 12 |
| products_reverse | Float32 | 0 | 31 |
| site_math | Float32 | 19 | 66 |
| force_virial_reduction | Float32 | 0 | 13 |
| sphericart_normalize | Float32 | 40 | 13 |
| sphericart_kernel | Float32 | 135 | 55 |
| sphericart_rescale | Float32 | 21 | 5 |

Imported cubic Horner contains packed `vmulpd`/`vaddpd` (or `ps`) arithmetic;
transform/envelope, binary interval selection and derivative setup remain scalar
per lane. The sampled Float64 spline uses scalar coefficient loads followed by
`vinsertf128` lane assembly, not hardware gathers. Custom analytic radial calls
remain lane-wise. Their exact recurrence is preserved rather than replaced by
an approximate spline.

SpheriCart uses contiguous temporary component arrays and compiler-vectorized
loops: the actual normalization, harmonic derivative and rescaling kernels
contain packed arithmetic in addition to scalar remainder work. The explicit
angular block also contains packed arithmetic, but its native code includes
allocation calls (`jl_alloc_genericmemory_unchecked`, `ijl_gc_small_alloc`),
consistent with measured host bytes. Steady prepared batched SpheriCart cases
report zero host allocation in these measurements.

Density accumulation and ragged product/reverse loops show predominantly scalar
arithmetic. Column-major R/Y keep neighbor samples contiguous, but each density
sum carries a recurrence and products follow variable index lists. The complete
site math contains both scalar and packed operations; it does not prove the
whole edge contraction is vectorized. The inspected sorted-force/virial scatter
remains scalar, with indexed writes protected by the existing colored schedule.
Geometry normalization does show compiler vectorization in the batched
SpheriCart path. Three-component `SVector` arithmetic and indirect scatter are
different layouts from the explicit neighbor SIMD lanes.

The CPU equivariant feature evaluator is a separate functional serial path. Source
inspection shows per-site radial/angular arrays and coupling/readout `mul!` calls;
this audit does not time or claim SIMD coverage for that feature path. Its existing
numerical/AD regression suite remains covered by the full package tests.

Measured priorities for later CPU work are density/edge contraction and ragged
product reverse; the public PI expansion makes product work substantial. Explore
SIMD across independent central sites only with compatible force scheduling, and
compare bounded preallocated angular temporaries before revisiting block fusion.
Retain batched SpheriCart as the measured baseline. No CPU kernel optimization is
introduced by this audit.

Raw reports and gzip-compressed LLVM/native text:
[dated artifacts](../benchmark/gpu_ace/results/2026-11-09-julia1.13.3-5eef1c1-dirty/).
Reproduce with `benchmark/gpu_ace/audit.jl` as documented in the benchmark README.
