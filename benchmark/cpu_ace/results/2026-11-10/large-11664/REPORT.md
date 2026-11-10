# CPU ACE stage-one measured results

Warmed BenchmarkTools medians. Ratios compare identical numerical models and boundaries; raw samples and allocation observations remain in the TOML files.

| Threads | Model | Atoms | Precision / width | Reference ms | Optimized ms | Incremental ratio | Historical ratio | Public ms | External public ms |
|---:|---|---:|---|---:|---:|---:|---:|---:|---:|
| 1 | medium-ti | 11664 | Float64 / 4 | 205.503 | 148.684 | 1.38x | unmeasured | 338.114 | 254.337 |
| 1 | large-tial | 11664 | Float64 / 4 | 6671.578 | 1604.239 | 4.16x | unmeasured | 1798.852 | 1383.335 |

External times are the stored oracle thread count/boundary, not the PFW thread count. Modern references are combined Float64 public E/F/V; Float32 rows show a Float64 sensitivity comparison. Legacy references are explicitly three separate analytic ACE1 calls with JuLIP native neighbor caching at fixed geometry. Do not interpret the prepared/external ratio as a pure kernel speedup.

| Model / precision | Algorithm | First-site product µs | Density µs | Radial µs | Angular µs | Force back-contraction µs |
|---|---|---:|---:|---:|---:|---:|

Isolated site stages perturb state/cache behavior and cannot be added to reconstruct total execution time. Preparation observations retain species packing and colored scheduling separately.

Under-sampled trials:

- `current-t1.toml` large-tial / 11664 / Float64 / w4 / reference / repeated_frames: 5 samples.

Omissions:

- `current-t1.toml` PACE: unsupported: Executable probed; equivalent exported spline interpolation and pair tables not validated in stage 1. No PACE comparison launched.

Source identities, package versions, scratch payloads, geometry/model hashes, absolute E/F/V errors and raw observations are in each TOML.
