# Recorded GPU ACE milestone 1 evidence

2026-11-09; base commit 5eef1c1e128385ce74950984256ca4590c0fce47, dirty
milestone working tree. CPU: Ryzen 9 9950X (`znver5`). GPU: Radeon RX 6800 XT,
gfx1030, 72 CUs. Julia 1.13.3, four process threads / one BLAS thread unless
specified; actual device/runtime/package versions are in the reports. The
Julia 1.10.12 CPU suite is recorded separately. Device nodes were hidden in the
restricted sandbox; these AMD results used the actual host-accessible GPU.

All timings are warmed, synchronized where applicable, `evals=1`, and collected
sequentially. Quick requests three samples / 0.2 seconds; full requests three /
0.1 seconds. The time budget can stop slow cases after one sample. These brief
runs demonstrate behavior and bottlenecks; they are not confidence intervals.
Throughput, individual samples, spread, host allocations and per-case errors are
retained. The quick workers' elapsed times include checks/compilation and are
about 18.6 s CPU / 62.2 s AMD; driver/probe startup is additional.

## Quick Float64 measurements

| Model | Atoms | CPU prepared (ms) | AMD staged prepared (ms) | CPU complete (ms) | AMD complete (ms) | Original oracle (ms) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| tial | 128 | 1.766 | 9.948 | 3.525 | 10.416 | 1.411 |
| tial | 11,664 | 50.169 | 29.442 | 231.730 | 32.940 | — |
| ti | 128 | 79.514 | 267.533 | 83.231 | 268.557 | 13.520 |
| tialv | 128 | 221.494 | 2475.375 | 228.683 | 2474.991 | 14.568 |

CPU uses the block angular path, four executor slots and colored scheduling.
Prepared GPU includes fresh sorted cells/connectivity at zero skin, readout
refresh, numerical E/F/V and status readback; model/inputs/scratch are resident.
Complete calls include boundary extraction, input update/upload, preparation,
execution and owning unitful host results. Initial explicit-strategy reports
(`quick`/`full`/`compact-scale`/`scale1024-*`) measured **the default staged public
calculator for every complete call**, including compact/baseline jobs. Their
strategy field labels prepared execution and trajectory. They predate automatic
compact selection. Later `auto-final`/`float32-final` reports explicitly warm the
public calculator with the same requested mapping, record `complete_strategy`,
actual/requested strategy, edge counts and full resident payload categories.
This reporting/preparation refinement changes no kernel mathematics.

Original TiAl is ACEpotentials 0.10.2 / Julia 1.13.3 and uses its combined public
E/F/V call. Ti and April08 Ti–Al–V are ACEpotentials 0.6.6 / Julia 1.10.12,
and sum the timings of three separate resident JuLIP property calls in one timed
closure. All original oracle producers use four Julia threads / one BLAS thread.
Those boundaries differ from prepared PFW execution; do not turn this table into
an unqualified library speedup claim. No oracle is timed at 11,664 atoms.

## Mapping experiment, 128 atoms, Float64

| Model | Staged (ms) | One-thread baseline (ms) | Compact (ms) |
| --- | ---: | ---: | ---: |
| tial | 10.205 | 13.928 | 1.428 |
| ti | 267.598 | 462.563 | 45.144 |
| tialv | 2484.901 | 4215.083 | 251.562 |

Compact has a strong small-system benefit, while the one-thread baseline loses.
For TiAl at 1,024 atoms, staged/compact medians are 11.066/3.845 ms (five samples,
0.5-second budget). At 11,664 atoms they are 29.442/29.052 ms in separate runs;
that small difference is not a robust compact advantage. The resulting automatic
rule is deliberately limited to RX 6800 XT, Float64, spline-only lmax ≤ 2 models
with ≤27 densities / ≤167 products and 128–1,024 atoms, selected at executor
construction. Others use staged by default and retain manual overrides. Larger
public models are calibrated only at 128 atoms here, so their compact benefit is
reported without extending the automatic policy to them.

Public-model force differences against independent Float64 originals are below
4e-14 eV/Å in the recorded Float64 fixtures. Float32 AMD differences reach about
1.5e-5 eV/Å, 0.002 eV virial and 0.083 eV total energy on these examples;
large fixed reference offsets contribute to total-energy rounding. GPU-versus-
rounded-CPU errors are separately recorded, not conflated with application
precision sensitivity. No universal Float32 suitability threshold is assigned.

The full run covers all three models at 128 atoms, Float32/Float64, lattice
parameters 3.4/3.2 Å, separate CPU one-/four-thread processes, all GPU mappings,
and deterministic five-frame cached-neighbor trajectories with forced rebuilds.
TiAl records two rebuilds, public models three, in the staged sequence. Larger
public-model full defaults (and 27,648-atom cases) were **not measured**; explicit
`--sizes 4` bounded the full matrix. The quick preset supplies the 10k+ workload.

CPU SIMD findings and limitations: [audit report](../../../../docs/gpu_ace_cpu_simd_audit.md).
`full/audit-*` preserves the original default-argument-wrapper inspection.
`audit-corrected/*` contains the corrected actual radial/site/SpheriCart kernels.
LLVM/native text is gzip-compressed without content edits. Static occurrences
must not be interpreted as runtime SIMD percentages.

## Reproduction

Run from the repository root. The environment snapshots rebase only the package's
local path to the repository; registry versions remain pinned. Use the matching
Julia major/minor for each snapshot, and `Pkg.instantiate()` if needed. Vendor
artifacts require the normal installed ROCm/device setup.

```sh
results=benchmark/gpu_ace/results/2026-11-09-julia1.13.3-5eef1c1-dirty
julia --project="$results/environments/cpu-1.13" -e 'using Pkg; Pkg.instantiate()'
julia --project="$results/environments/amd-1.13" -e 'using Pkg; Pkg.instantiate()'
julia benchmark/gpu_ace/suite.jl --preset quick --backend both --data /tmp/pfw-public-ace --cpu-project "$results/environments/cpu-1.13" --amd-project "$results/environments/amd-1.13" --samples 3 --seconds 0.2 --output /tmp/reproduced-quick
julia benchmark/gpu_ace/suite.jl --preset full --backend both --data /tmp/pfw-public-ace --cpu-project "$results/environments/cpu-1.13" --amd-project "$results/environments/amd-1.13" --samples 3 --seconds 0.1 --sizes 4 --output /tmp/reproduced-full
JULIA_LOAD_PATH="@:$PWD:@stdlib" julia --project="$results/environments/amd-1.13" -t4 test/gpu_ace.jl --amdgpu
JULIA_LOAD_PATH="@:$PWD:@stdlib" julia --project="$results/environments/cpu-1.13" -t4 test/runtests.jl
```

First acquire/export/verify public models using the benchmark README. The saved
source JSONs and large generated exports remain in the caller's cache; checksums,
original oracle structures/predictions, raw archive metadata/notices and
verification hashes are in `validation`. Dataset attribution and licensing are
in `../../THIRD_PARTY.md`, separate from package code licensing. Reports retain
original absolute run/log paths in `/tmp`; their copied siblings here contain
the corresponding artifacts. See `../../REPORTS.md` for field definitions.

GPU AD/training/features, skew/multiple-image geometry and other vendors remain
outside validation. Device allocation peaks/register spills are unavailable;
array payload counts and host allocations are different measurements. Large
legacy analytic radials repeat channel work in the force pass and remain a
performance limitation; no speed requirement was substituted for correctness.

Final checks after adding dispatch/report refinements: actual RX ACE 758 assertions,
portable kernels 746 assertions; automatic TiAl runs at 128/1,024/11,664 atoms;
Float32 implementation energy/force/virial checks for all three models; and a
CPU Float32 worker invoked directly through the rebased environment snapshot.
The final CPU worker is `cpu-final.toml`, with its own one-thread boundary.
