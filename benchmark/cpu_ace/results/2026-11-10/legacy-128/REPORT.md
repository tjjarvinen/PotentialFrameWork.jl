# CPU ACE stage-one measured results

Warmed BenchmarkTools medians. Ratios compare identical numerical models and boundaries; raw samples and allocation observations remain in the TOML files.

| Threads | Model | Atoms | Precision / width | Reference ms | Optimized ms | Incremental ratio | Historical ratio | Public ms | External public ms |
|---:|---|---:|---|---:|---:|---:|---:|---:|---:|
| 1 | ti | 128 | Float64 / 4 | 74.664 | 24.805 | 3.01x | 3.00x | 28.238 | 47.287 |
| 1 | ti | 128 | Float32 / 4 | 73.463 | 23.686 | 3.10x | 3.16x | 24.187 | 47.287 |
| 1 | tialv | 128 | Float64 / 4 | 217.118 | 60.118 | 3.61x | 3.61x | 64.925 | 51.759 |
| 1 | tialv | 128 | Float32 / 4 | 214.236 | 57.837 | 3.70x | 3.80x | 59.732 | 51.759 |

External times are the stored oracle thread count/boundary, not the PFW thread count. Modern references are combined Float64 public E/F/V; Float32 rows show a Float64 sensitivity comparison. Legacy references are explicitly three separate analytic ACE1 calls with JuLIP native neighbor caching at fixed geometry. Do not interpret the prepared/external ratio as a pure kernel speedup.

| Model / precision | Algorithm | First-site product µs | Density µs | Radial µs | Angular µs | Force back-contraction µs |
|---|---|---:|---:|---:|---:|---:|
| ti / Float64 / w4 | reference | 526.573 | 33.135 | 10.900 | 3.070 | 22.120 |
| ti / Float64 / w4 | optimized | 105.141 | 32.635 | 11.000 | 3.050 | 21.765 |
| ti / Float32 / w4 | reference | 487.987 | 33.210 | 10.310 | 1.765 | 20.055 |
| ti / Float32 / w4 | optimized | 104.500 | 33.000 | 10.235 | 1.765 | 19.855 |
| tialv / Float64 / w4 | reference | 1734.164 | 42.116 | 12.316 | 1.605 | 30.555 |
| tialv / Float64 / w4 | optimized | 370.247 | 41.760 | 12.345 | 1.585 | 29.010 |
| tialv / Float32 / w4 | reference | 1569.293 | 43.306 | 11.281 | 0.865 | 27.365 |
| tialv / Float32 / w4 | optimized | 348.637 | 42.100 | 11.395 | 0.865 | 27.965 |

Isolated site stages perturb state/cache behavior and cannot be added to reconstruct total execution time. Preparation observations retain species packing and colored scheduling separately.

Under-sampled trials:

None for the three reported whole-call boundaries.

Omissions:

- `current-t1.toml` PACE: unsupported: Executable probed; equivalent exported spline interpolation and pair tables not validated in stage 1. No PACE comparison launched.
- `historical-t1.toml` PACE: unsupported: Executable probed; equivalent exported spline interpolation and pair tables not validated in stage 1. No PACE comparison launched.

Source identities, package versions, scratch payloads, geometry/model hashes, absolute E/F/V errors and raw observations are in each TOML.
