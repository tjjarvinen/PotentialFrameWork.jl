# CPU ACE stage-one measured results

Warmed BenchmarkTools medians. Ratios compare identical numerical models and boundaries; raw samples and allocation observations remain in the TOML files.

| Threads | Model | Atoms | Precision / width | Reference ms | Optimized ms | Incremental ratio | Historical ratio | Public ms | External public ms |
|---:|---|---:|---|---:|---:|---:|---:|---:|---:|
| 1 | small-tial | 128 | Float64 / 4 | 1.743 | 1.317 | 1.32x | 1.33x | 3.043 | 2.362 |
| 1 | small-tial | 128 | Float32 / 4 | 1.649 | 1.225 | 1.35x | 1.35x | 1.533 | 2.362 |
| 1 | small-tial | 11664 | Float64 / 4 | 159.648 | 120.121 | 1.33x | 1.32x | 317.670 | 232.715 |
| 1 | small-tial | 11664 | Float32 / 4 | 150.148 | 112.221 | 1.34x | 1.33x | 145.500 | 232.715 |
| 1 | medium-ti | 128 | Float64 / 4 | 2.239 | 1.621 | 1.38x | 1.38x | 3.134 | 2.551 |
| 1 | medium-ti | 128 | Float32 / 4 | 2.131 | 1.488 | 1.43x | 1.40x | 1.760 | 2.551 |
| 1 | large-tial | 128 | Float64 / 4 | 72.907 | 17.584 | 4.15x | 4.15x | 19.392 | 16.638 |
| 1 | large-tial | 128 | Float32 / 4 | 71.170 | 17.345 | 4.10x | 4.21x | 17.744 | 16.638 |

External times are the stored oracle thread count/boundary, not the PFW thread count. Modern references are combined Float64 public E/F/V; Float32 rows show a Float64 sensitivity comparison. Legacy references are explicitly three separate analytic ACE1 calls with JuLIP native neighbor caching at fixed geometry. Do not interpret the prepared/external ratio as a pure kernel speedup.

| Model / precision | Algorithm | First-site product µs | Density µs | Radial µs | Angular µs | Force back-contraction µs |
|---|---|---:|---:|---:|---:|---:|
| small-tial / Float64 / w4 | reference | 4.180 | 1.010 | 6.135 | 0.390 | 1.070 |
| small-tial / Float64 / w4 | optimized | 0.930 | 1.000 | 6.000 | 0.390 | 1.050 |
| small-tial / Float32 / w4 | reference | 4.085 | 1.005 | 5.515 | 0.285 | 0.985 |
| small-tial / Float32 / w4 | optimized | 0.940 | 0.980 | 5.530 | 0.270 | 0.980 |
| medium-ti / Float64 / w4 | reference | 5.890 | 1.660 | 6.270 | 0.850 | 1.470 |
| medium-ti / Float64 / w4 | optimized | 1.395 | 1.650 | 6.335 | 0.825 | 1.490 |
| medium-ti / Float32 / w4 | reference | 5.815 | 1.660 | 5.890 | 0.585 | 1.375 |
| medium-ti / Float32 / w4 | optimized | 1.315 | 1.650 | 5.860 | 0.585 | 1.370 |
| large-tial / Float64 / w4 | reference | 530.722 | 5.910 | 11.445 | 1.605 | 5.320 |
| large-tial / Float64 / w4 | optimized | 98.120 | 5.895 | 11.330 | 1.550 | 5.205 |
| large-tial / Float32 / w4 | reference | 578.328 | 5.880 | 10.715 | 0.975 | 4.950 |
| large-tial / Float32 / w4 | optimized | 100.435 | 5.820 | 10.765 | 0.960 | 4.935 |

Isolated site stages perturb state/cache behavior and cannot be added to reconstruct total execution time. Preparation observations retain species packing and colored scheduling separately.

Under-sampled trials:

None for the three reported whole-call boundaries.

Omissions:

- `current-t1.toml` PACE: unsupported: Executable probed; equivalent exported spline interpolation and pair tables not validated in stage 1. No PACE comparison launched.
- `historical-t1.toml` PACE: unsupported: Executable probed; equivalent exported spline interpolation and pair tables not validated in stage 1. No PACE comparison launched.

Source identities, package versions, scratch payloads, geometry/model hashes, absolute E/F/V errors and raw observations are in each TOML.
