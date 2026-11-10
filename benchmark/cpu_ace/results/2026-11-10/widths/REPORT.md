# CPU ACE stage-one measured results

Warmed BenchmarkTools medians. Ratios compare identical numerical models and boundaries; raw samples and allocation observations remain in the TOML files.

| Threads | Model | Atoms | Precision / width | Reference ms | Optimized ms | Incremental ratio | Historical ratio | Public ms | External public ms |
|---:|---|---:|---|---:|---:|---:|---:|---:|---:|
| 1 | small-tial | 128 | Float64 / 1 | 1.647 | 1.223 | 1.35x | unmeasured | 2.909 | 2.362 |
| 1 | small-tial | 128 | Float64 / 8 | 2.261 | 1.818 | 1.24x | unmeasured | 3.517 | 2.362 |
| 1 | small-tial | 128 | Float32 / 1 | 1.468 | 1.059 | 1.39x | unmeasured | 1.322 | 2.362 |
| 1 | small-tial | 128 | Float32 / 8 | 2.036 | 1.622 | 1.26x | unmeasured | 1.911 | 2.362 |

External times are the stored oracle thread count/boundary, not the PFW thread count. Modern references are combined Float64 public E/F/V; Float32 rows show a Float64 sensitivity comparison. Legacy references are explicitly three separate analytic ACE1 calls with JuLIP native neighbor caching at fixed geometry. Do not interpret the prepared/external ratio as a pure kernel speedup.

| Model / precision | Algorithm | First-site product µs | Density µs | Radial µs | Angular µs | Force back-contraction µs |
|---|---|---:|---:|---:|---:|---:|
| small-tial / Float64 / w1 | reference | 4.340 | 0.920 | 7.300 | 0.355 | 1.055 |
| small-tial / Float64 / w1 | optimized | 0.955 | 0.900 | 7.420 | 0.350 | 1.055 |
| small-tial / Float64 / w8 | reference | 4.325 | 1.150 | 9.205 | 0.370 | 1.065 |
| small-tial / Float64 / w8 | optimized | 0.960 | 1.150 | 9.225 | 0.365 | 1.050 |
| small-tial / Float32 / w1 | reference | 4.080 | 0.920 | 6.610 | 0.315 | 0.970 |
| small-tial / Float32 / w1 | optimized | 0.960 | 0.900 | 6.540 | 0.310 | 0.970 |
| small-tial / Float32 / w8 | reference | 4.320 | 1.150 | 7.885 | 0.250 | 0.985 |
| small-tial / Float32 / w8 | optimized | 0.960 | 1.150 | 7.980 | 0.250 | 0.980 |

Isolated site stages perturb state/cache behavior and cannot be added to reconstruct total execution time. Preparation observations retain species packing and colored scheduling separately.

Under-sampled trials:

None for the three reported whole-call boundaries.

Omissions:

- `current-t1.toml` PACE: unsupported: Executable probed; equivalent exported spline interpolation and pair tables not validated in stage 1. No PACE comparison launched.

Source identities, package versions, scratch payloads, geometry/model hashes, absolute E/F/V errors and raw observations are in each TOML.
