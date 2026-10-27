# GPU phase 1 results — 27 October 2026

Base commit: `a6b3dd56e2560fe1339f0f41464da052f5cb7ff7`.
The user's appended AGENTS.md, GPU phase documents and local ACE reference data
were preserved. Source, examples and tests are local changes, not a remote reset.

The real-device suite passed **1,539 assertions** on the AMD Radeon RX 6800 XT
(gfx1030, 15.984 GiB, reported wavefront 32), using Julia 1.13.3. This includes
Float32/64 and precision overrides, all requested property subsets, independent
O(N²) references on converted inputs, finite difference force/strain, strict
named/isotopic identity through actual device accessors, updated parameters,
owning outputs, sparse/open/mixed/periodic cells, wrapped cell alias removal,
large/tail populations, alternate launches, singular origins in the convenience
and device-numerical paths, skin threshold/approaching-pair/invalidation cases.
Both strategies pass. No GPU AD or GPU ACE/MACE is claimed.

Hardware/software diagnostics and the real-device test summaries are saved in
[GPU-PAIR-HARDWARE.txt](GPU-PAIR-HARDWARE.txt). The host is Linux
7.2.11-200.fc44.x86_64 on an AMD Ryzen 9 9950X (16 cores/32 logical CPUs), with
HIP 7.1.52802. GPU validation used AMDGPU 2.10.0, GPUCompiler 2.18.0, GPUArrays
12.1.1, LLVM.jl 10.4.2, KA 0.10.3 and AK 0.5.10. Other relevant versions are
AtomsBase 0.5.3, AtomsCalculators 0.2.3, AtomsBuilder 0.2.3, AtomsSystems 0.1.2,
SIMD 3.8.1, StaticArrays 1.9.24 and Unitful 1.29.0. All application dependency
versions are in the measurement TOML. Benchmarks used 32 Julia threads; full CPU
regressions use four. The integrated GPU was enumerated but not used.

## Fresh-cell timings

Warm medians in milliseconds, seven synchronized repetitions. The numerical
column excludes input upload and host force output, but includes fresh cells,
record sorting, search, callback mathematics, reduction, a scalar validity read
and owning device result copies. The public column includes AtomsSystems
extraction/upload and owning unitful results. The CPU column uses the existing
parallel width-4 public CPU evaluator on the same converted geometry/parameters.
No validation process ran concurrently with these final timing samples.

| Atoms | Strategy | Float32 argon E/F/W numerical | Public GPU | Public CPU |
| ---: | --- | ---: | ---: | ---: |
| 10,976 | central | 0.479 | 0.728 | 4.741 |
| 10,976 | candidate | 1.395 | 1.639 | 4.595 |
| 97,556 | central | 1.477 | 2.906 | 28.655 |
| 97,556 | candidate | 11.693 | 13.086 | 31.667 |

The experimental candidate strategy is slower in these fresh-cell samples and is
retained as requested. This does not establish a general speed ratio across
models, systems or hardware. The complete 48-case dataset includes energy-only
and full E/F/W, Float32/64, argon LJ, a mixed Ar/Ne SmoothPair and an open-z slab.
It records first preparation/execution (compilation/capacity growth when needed),
min/max, stage profiles, useful atom-step rate, host allocations, device buffers
and ownership-copy sizes. The source fixes seeds, perturbations, geometry and
species fractions. Stage profiles synchronize separately and are not additive
estimates of the asynchronous total. The keys-and-sort profile includes key
regeneration so every sorting sample receives the same unsorted records.

LJ is shifted at 8.5 Å with epsilon=0.0103 eV and sigma=3.4 Å; FCC lattice
constant is 5.26 Å, density 0.0274854441 atoms/Å³, bounded rattle 0.03 Å.
Replication 14³/29³ gives 10,976/97,556 atoms, cell lengths 73.64/152.54 Å.
The alternating mixed/slab species fractions are 0.5/0.5. Default launch is 64
lanes and four centers per tile. All three requested outputs pass reference
gates before timing, including strain-convention virials.

## Numerical error and storage

Focused checks use absolute/relative tolerances 3e-5 for Float32 and 3e-11 for
Float64; these accommodate input rounding and different reduction order for the
tested SmoothPair/LJ scales. Independent sums use the explicitly converted input
coordinates and parameter precision, with double arithmetic, separating input
rounding from execution error. Large benchmark gates use 2e-4/2e-10 respectively,
including a separate absolute maximum-force gate. Measured maxima over all
benchmark cases were:

| Compute | Max per-atom force error vs same-precision CPU (eV/Å) | Max per-atom force error vs original Float64 input (eV/Å) | Max energy error/atom vs original (eV) |
| --- | ---: | ---: | ---: |
| Float32 | 2.22e-06 | 5.59e-06 | 6.4e-09 |
| Float64 | 4.08e-15 | 4.08e-15 | 1.55e-17 |

Per-atom errors are Euclidean vector norms, recorded through the deterministic
`benchmark/gpu_pair_force_errors.jl` replay without changing timing fields. The
main benchmark records these fields directly on subsequent runs. Maximum
component errors are also retained in the dataset.

At 97,556 atoms the Float32 argon full E/F/W case holds 22,643,320 bytes of counted main/upload/algorithm-scratch/owning-result device buffers. No full
N-atom force buffer per workgroup or dense distance matrix is used. Numerical
results deliberately allocate owning outputs. Host allocation counts describe
Julia wrappers, not GPU allocator transactions; no allocation-free claim is made.

## Evolving skin trajectory

Forty nonuniform displacement steps (component amplitude 0.005 Å/step) start from
the same 10,976-atom perturbed argon structure. Both rebuild and reuse specializations
are warmed before timing. Input upload and fresh-reference comparisons are outside
these numerical timings; every step is checked against fresh GPU execution.

| Strategy | Skin (Å) | Rebuild steps / 40 | Mean ms/step | Median reuse ms |
| --- | ---: | ---: | ---: | ---: |
| central | 0.0 | 40 | 0.455 | — |
| central | 0.5 | 1 | 0.256 | 0.169 |
| candidate | 0.0 | 40 | 1.371 | — |
| candidate | 0.5 | 1 | 0.242 | 0.220 |

Skin=0.5 uses 7,304,328 additional connectivity/reference/scan/reduction bytes in
this trajectory. Rebuild steps are substantially more expensive than reused-list
steps; all individual samples are in `benchmark/gpu_skin_*.toml`. Skin improves
this particular trajectory; it is not automatically selected or tuned.

## Compatibility and reproduction

Full CPU/interface/AD results and final counts are recorded in STATUS.md.
Standalone core import passed with neither AMDGPU nor CUDA loaded. KA 0.9/AK 0.4
retain the existing SpheriCart 0.2 environment, while the real-device application
uses KA 0.10/AK 0.5. Both KA CPU backend generations were checked; POCL's initial
struct reinterpretation failure was fixed with explicit species field decoding.
Earlier AMDGPU 2.8 attempts failed on Julia 1.12/1.13 with sort-kernel InvalidIR,
and a Julia 1.10 LLVM compilation attempt crashed. These combinations are not
claimed as supported GPU configurations. Core Julia compatibility stays 1.10.10;
1.10.12 is the tested installed patch release.

See [../gpu_pairs.md](../gpu_pairs.md) for commands, callback certification,
position/parameter refresh and input/output lifetimes. GPU geometry is restricted
to single-image orthorhombic/open/mixed cells with strict cutoff+skin eligibility;
skew and multi-image GPU requests reject. CPU skew/multi-image and functional AD
remain supported. Other vendors, GPU AD, GPU ACE/MACE and optimized Newton/atomic
reductions remain outside this delivery. The G1 radial correctness gate uses production device-built sorted cells and
ranges. No transferred-CPU-connectivity development path is included.
