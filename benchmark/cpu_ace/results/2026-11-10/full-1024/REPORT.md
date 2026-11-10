# CPU ACE stage-one measured results

Warmed BenchmarkTools medians. Ratios compare identical numerical models and boundaries; raw samples and allocation observations remain in the TOML files.

| Threads | Model | Atoms | Precision / width | Reference ms | Optimized ms | Incremental ratio | Historical ratio | Public ms | External public ms |
|---:|---|---:|---|---:|---:|---:|---:|---:|---:|
| 1 | small-tial | 1024 | Float64 / 4 | 14.041 | 10.607 | 1.32x | 1.33x | 25.671 | 19.226 |
| 1 | small-tial | 1024 | Float32 / 4 | 13.092 | 9.792 | 1.34x | 1.36x | 12.559 | 19.226 |
| 1 | medium-ti | 1024 | Float64 / 4 | 17.847 | 13.014 | 1.37x | 1.38x | 26.366 | 20.908 |
| 1 | medium-ti | 1024 | Float32 / 4 | 16.955 | 12.199 | 1.39x | 1.39x | 14.692 | 20.908 |
| 1 | large-tial | 1024 | Float64 / 4 | 578.571 | 137.960 | 4.19x | 4.27x | 154.870 | 122.978 |
| 1 | large-tial | 1024 | Float32 / 4 | 570.991 | 136.133 | 4.19x | 4.25x | 139.246 | 122.978 |
| 16 | small-tial | 1024 | Float64 / 4 | 3.741 | 2.698 | 1.39x | 1.32x | 18.581 | 19.226 |
| 16 | small-tial | 1024 | Float32 / 4 | 2.907 | 2.446 | 1.19x | 1.22x | 5.516 | 19.226 |
| 16 | medium-ti | 1024 | Float64 / 4 | 4.217 | 3.197 | 1.32x | 1.32x | 17.410 | 20.908 |
| 16 | medium-ti | 1024 | Float32 / 4 | 3.268 | 2.591 | 1.26x | 1.34x | 5.984 | 20.908 |
| 16 | large-tial | 1024 | Float64 / 4 | 133.641 | 26.666 | 5.01x | 5.97x | 43.732 | 122.978 |
| 16 | large-tial | 1024 | Float32 / 4 | 89.318 | 25.241 | 3.54x | 3.58x | 29.694 | 122.978 |
| 4 | small-tial | 1024 | Float64 / 4 | 4.278 | 3.345 | 1.28x | 1.28x | 18.991 | 19.226 |
| 4 | small-tial | 1024 | Float32 / 4 | 3.965 | 3.034 | 1.31x | 1.30x | 6.142 | 19.226 |
| 4 | medium-ti | 1024 | Float64 / 4 | 5.346 | 3.945 | 1.36x | 1.43x | 18.040 | 20.908 |
| 4 | medium-ti | 1024 | Float32 / 4 | 4.914 | 3.578 | 1.37x | 1.43x | 6.923 | 20.908 |
| 4 | large-tial | 1024 | Float64 / 4 | 166.606 | 39.572 | 4.21x | 6.54x | 55.509 | 122.978 |
| 4 | large-tial | 1024 | Float32 / 4 | 164.453 | 39.524 | 4.16x | 4.17x | 43.158 | 122.978 |

External times are the stored oracle thread count/boundary, not the PFW thread count. Modern references are combined Float64 public E/F/V; Float32 rows show a Float64 sensitivity comparison. Legacy references are explicitly three separate analytic ACE1 calls with JuLIP native neighbor caching at fixed geometry. Do not interpret the prepared/external ratio as a pure kernel speedup.

| Model / precision | Algorithm | First-site product µs | Density µs | Radial µs | Angular µs | Force back-contraction µs |
|---|---|---:|---:|---:|---:|---:|
| small-tial / Float64 / w4 | reference | 4.160 | 1.000 | 5.785 | 0.390 | 1.065 |
| small-tial / Float64 / w4 | optimized | 0.950 | 0.980 | 5.885 | 0.380 | 1.050 |
| small-tial / Float32 / w4 | reference | 4.075 | 1.010 | 5.425 | 0.295 | 0.990 |
| small-tial / Float32 / w4 | optimized | 0.950 | 0.990 | 5.540 | 0.290 | 0.980 |
| medium-ti / Float64 / w4 | reference | 5.895 | 1.660 | 6.140 | 0.830 | 1.480 |
| medium-ti / Float64 / w4 | optimized | 1.320 | 1.650 | 6.185 | 0.830 | 1.490 |
| medium-ti / Float32 / w4 | reference | 5.831 | 1.695 | 5.805 | 0.585 | 1.380 |
| medium-ti / Float32 / w4 | optimized | 1.390 | 1.650 | 5.965 | 0.565 | 1.430 |
| large-tial / Float64 / w4 | reference | 533.528 | 5.900 | 11.410 | 1.590 | 5.225 |
| large-tial / Float64 / w4 | optimized | 98.535 | 5.825 | 11.025 | 1.535 | 5.140 |
| large-tial / Float32 / w4 | reference | 525.903 | 5.950 | 10.646 | 0.965 | 4.750 |
| large-tial / Float32 / w4 | optimized | 98.271 | 5.890 | 10.720 | 0.940 | 4.950 |

Isolated site stages perturb state/cache behavior and cannot be added to reconstruct total execution time. Preparation observations retain species packing and colored scheduling separately.

Under-sampled trials:

None for the three reported whole-call boundaries.

Omissions:

- `current-t1.toml` PACE: unsupported: Executable probed; equivalent exported spline interpolation and pair tables not validated in stage 1. No PACE comparison launched.
- `current-t16.toml` PACE: unsupported: Executable probed; equivalent exported spline interpolation and pair tables not validated in stage 1. No PACE comparison launched.
- `current-t4.toml` PACE: unsupported: Executable probed; equivalent exported spline interpolation and pair tables not validated in stage 1. No PACE comparison launched.
- `historical-t1.toml` PACE: unsupported: Executable probed; equivalent exported spline interpolation and pair tables not validated in stage 1. No PACE comparison launched.
- `historical-t16.toml` PACE: unsupported: Executable probed; equivalent exported spline interpolation and pair tables not validated in stage 1. No PACE comparison launched.
- `historical-t4.toml` PACE: unsupported: Executable probed; equivalent exported spline interpolation and pair tables not validated in stage 1. No PACE comparison launched.

Source identities, package versions, scratch payloads, geometry/model hashes, absolute E/F/V errors and raw observations are in each TOML.
