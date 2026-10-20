# Recorded Argon Lennard-Jones comparison

Recorded environment: `Julia=1.13.2 default_threads=1 CPU=AMD Ryzen 9 9950X 16-Core Processor compiler_target=znver5 affinity=uncontrolled precision=Float64 warmups=3 samples=7 cells=[4, 6, 10, 16]`

Full calls use `AtomsCalculators.energy_forces_virial` with the matched inverse-distance cutoff tilt.
Reference: EmpiricalPotentials 0.2.6; systems: AtomsSystems 0.1.2 perturbed periodic FCC Ar.
Worker processes ran sequentially, with 3 warmups and 7 measured samples per case.
CPU affinity was uncontrolled. Exact versions/hardware are in the metadata files.
All values below are milliseconds; brackets show minimum/maximum around the median.
`EP/PFW4` compares implementations; `W4/W8` compares vector widths (>1 favours W=8).
Ratios use measured medians, not universal speedup claims or confidence intervals.

## Complete public calls, rebuilding topology each call

| N | Threads | EP threaded median [min, max] | PFW W=4 phased median [min, max] | EP/PFW4 | PFW W=8 phased median [min, max] | W4/W8 |
|---:|---:|---:|---:|---:|---:|---:|
| 256 | 1 | 2.798 [2.779, 2.855] | 0.221 [0.213, 0.347] | 12.64 | 0.209 [0.197, 0.337] | 1.059 |
| 256 | 2 | 2.709 [2.390, 2.819] | 0.257 [0.243, 0.334] | 10.55 | 0.230 [0.206, 0.297] | 1.118 |
| 256 | 4 | 2.299 [2.167, 2.349] | 0.284 [0.253, 0.434] | 8.09 | 0.240 [0.215, 0.290] | 1.184 |
| 256 | 8 | 1.973 [1.959, 2.068] | 0.275 [0.256, 0.387] | 7.18 | 0.258 [0.247, 0.398] | 1.063 |
| 864 | 1 | 9.783 [9.757, 9.803] | 1.431 [1.426, 1.496] | 6.84 | 1.545 [1.531, 1.753] | 0.926 |
| 864 | 2 | 8.217 [8.043, 8.790] | 1.463 [1.454, 1.843] | 5.62 | 1.542 [1.501, 1.690] | 0.949 |
| 864 | 4 | 7.550 [7.145, 8.110] | 1.481 [1.471, 1.532] | 5.10 | 1.509 [1.497, 1.564] | 0.981 |
| 864 | 8 | 7.013 [6.717, 8.270] | 1.516 [1.513, 1.584] | 4.63 | 1.543 [1.534, 1.705] | 0.982 |
| 4000 | 1 | 37.035 [35.061, 46.048] | 5.270 [5.262, 5.295] | 7.03 | 5.309 [5.305, 5.339] | 0.993 |
| 4000 | 2 | 19.434 [18.400, 25.245] | 3.400 [3.356, 3.890] | 5.72 | 3.417 [3.395, 3.967] | 0.995 |
| 4000 | 4 | 11.452 [9.762, 15.319] | 2.024 [2.006, 2.062] | 5.66 | 2.035 [2.024, 2.206] | 0.995 |
| 4000 | 8 | 10.445 [7.784, 33.279] | 1.551 [1.536, 1.880] | 6.74 | 1.687 [1.611, 1.861] | 0.919 |
| 16384 | 1 | 171.789 [159.884, 213.255] | 22.523 [22.401, 22.564] | 7.63 | 22.658 [22.620, 22.903] | 0.994 |
| 16384 | 2 | 93.844 [80.865, 123.577] | 13.538 [13.245, 13.972] | 6.93 | 13.421 [13.226, 14.386] | 1.009 |
| 16384 | 4 | 48.065 [45.428, 74.217] | 8.462 [8.384, 9.575] | 5.68 | 8.841 [8.495, 9.380] | 0.957 |
| 16384 | 8 | 38.132 [27.179, 63.380] | 5.749 [5.711, 6.255] | 6.63 | 6.219 [5.764, 7.098] | 0.924 |

## Unchanged-geometry diagnostics, setup excluded

| N | Threads | EP threaded median [min, max] | PFW W=4 phased median [min, max] | EP/PFW4 | PFW W=8 phased median [min, max] | W4/W8 |
|---:|---:|---:|---:|---:|---:|---:|
| 256 | 1 | 1.174 [1.155, 1.221] | 0.206 [0.199, 0.260] | 5.69 | 0.185 [0.182, 0.249] | 1.115 |
| 256 | 2 | 0.596 [0.590, 0.615] | 0.220 [0.210, 0.391] | 2.71 | 0.204 [0.194, 0.251] | 1.076 |
| 256 | 4 | 0.309 [0.305, 0.638] | 0.237 [0.213, 0.275] | 1.30 | 0.224 [0.199, 0.261] | 1.057 |
| 256 | 8 | 0.280 [0.272, 0.306] | 0.255 [0.232, 0.304] | 1.10 | 0.235 [0.224, 0.286] | 1.085 |
| 864 | 1 | 3.971 [3.958, 4.029] | 1.371 [1.369, 1.418] | 2.90 | 1.405 [1.400, 1.455] | 0.976 |
| 864 | 2 | 2.052 [2.027, 2.092] | 1.425 [1.415, 1.463] | 1.44 | 1.452 [1.427, 1.502] | 0.981 |
| 864 | 4 | 1.068 [1.035, 1.575] | 1.427 [1.417, 1.467] | 0.75 | 1.458 [1.454, 1.503] | 0.979 |
| 864 | 8 | 0.594 [0.571, 0.621] | 1.470 [1.463, 1.527] | 0.40 | 1.505 [1.482, 1.550] | 0.977 |
| 4000 | 1 | 20.865 [18.707, 23.146] | 4.889 [4.869, 4.911] | 4.27 | 4.949 [4.931, 4.974] | 0.988 |
| 4000 | 2 | 10.959 [9.470, 13.755] | 3.050 [3.006, 3.310] | 3.59 | 3.052 [3.029, 3.278] | 0.999 |
| 4000 | 4 | 5.765 [4.913, 9.593] | 1.672 [1.661, 2.043] | 3.45 | 1.687 [1.668, 1.826] | 0.991 |
| 4000 | 8 | 4.311 [4.124, 8.759] | 1.234 [1.210, 1.362] | 3.49 | 1.165 [1.153, 1.193] | 1.059 |
| 16384 | 1 | 82.608 [77.619, 88.450] | 21.021 [20.997, 21.054] | 3.93 | 21.139 [21.120, 21.177] | 0.994 |
| 16384 | 2 | 44.192 [40.979, 49.732] | 12.239 [12.100, 12.584] | 3.61 | 11.800 [11.786, 12.316] | 1.037 |
| 16384 | 4 | 24.013 [21.865, 34.553] | 7.178 [7.096, 7.316] | 3.35 | 7.276 [7.062, 9.783] | 0.987 |
| 16384 | 8 | 15.823 [14.565, 29.411] | 4.832 [4.529, 5.493] | 3.28 | 4.941 [4.297, 6.094] | 0.978 |

## Serial vector-width scaling

These rows use the 1-thread process; the executor is serial.
W=1 uses scalar pair arithmetic; W=4/8 use explicit SIMD vectors (Float64).
W1/W4 and W1/W8 compare against the scalar path; values >1 favour SIMD.

| N | Regime | W=1 median | W=4 median | W=8 median | W1/W4 | W1/W8 | W4/W8 |
|---:|---|---:|---:|---:|---:|---:|---:|
| 256 | rebuild | 0.185 | 0.214 | 0.196 | 0.864 | 0.945 | 1.095 |
| 256 | fixed | 0.166 | 0.201 | 0.187 | 0.829 | 0.890 | 1.075 |
| 864 | rebuild | 1.458 | 1.411 | 1.454 | 1.034 | 1.003 | 0.970 |
| 864 | fixed | 1.397 | 1.353 | 1.392 | 1.032 | 1.004 | 0.972 |
| 4000 | rebuild | 5.200 | 5.237 | 5.315 | 0.993 | 0.978 | 0.985 |
| 4000 | fixed | 4.831 | 4.869 | 4.969 | 0.992 | 0.972 | 0.980 |
| 16384 | rebuild | 22.284 | 22.186 | 22.524 | 1.004 | 0.989 | 0.985 |
| 16384 | fixed | 20.835 | 20.750 | 21.044 | 1.004 | 0.990 | 0.986 |

## Thread scaling at N=16384

Each ratio compares against the same implementation/regime at 1 thread(s).

| Threads | EP rebuild | PFW W=4 rebuild | PFW W=8 rebuild | EP fixed | PFW W=4 fixed | PFW W=8 fixed |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| 2 | 1.83 | 1.66 | 1.69 | 1.87 | 1.72 | 1.79 |
| 4 | 3.57 | 2.66 | 2.56 | 3.44 | 2.93 | 2.91 |
| 8 | 4.51 | 3.92 | 3.64 | 5.22 | 4.35 | 4.28 |

## Allocated bytes per complete call at N=16384

| Threads | EP threaded allocated bytes | PFW W=4 allocated bytes | PFW sorted force bytes | PFW phases |
|---:|---:|---:|---:|---:|
| 1 | 177,751,808 | 3,236,288 | 393,216 | 27 |
| 2 | 178,937,896 | 3,246,224 | 393,216 | 27 |
| 4 | 181,306,200 | 3,265,664 | 393,216 | 27 |
| 8 | 186,042,456 | 3,304,544 | 393,216 | 27 |

## Numerical agreement

All 224 measured cases passed energy, force and virial checks before/after timing.
Largest recorded absolute differences across cases:

- energy: 5.912e-12 eV
- force component: 4.025e-16 eV/Å
- virial component: 5.059e-12 eV

This compares complete implementations: a single-species analytic callback vs a general LJ implementation,
directed vs unordered pair ownership, automatic vs analytic
derivatives/precomputed cutoff coefficients, list construction, allocation, scheduling and reduction.
The fixed-topology table excludes setup, requires unchanged coordinates/cell/species, and is not
a public Verlet reuse claim. Small systems and phase occupancy can limit thread scaling.
Allocated bytes measure Julia allocations, not peak RSS; reference site assembly uses full force
partials per task, whereas PFW's sorted force workspace is independent of thread count.

See [methodology](README.md), [CSV records](results/), [time scaling](scaling.svg),
[vector-width scaling](vector-width-scaling.svg), [native-code inspection](codegen/README.md) and per-thread metadata.
