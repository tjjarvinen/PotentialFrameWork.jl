# GPU ACE milestone 1 workloads

Run commands from the repository root. These environments are intentionally
separate: SpheriCart 0.2.4 pins older KA, while the tested Julia 1.13 AMD stack
uses AMDGPU 2.10 / KA 0.10. The legacy oracle requires ACEpotentials 0.6.6,
ACEfit 0.1/0.2 and older AtomsBase; it cannot share the application environment.
No normal package import or test downloads public models.

```sh
julia --project=benchmark/gpu_ace -e 'using Pkg; Pkg.develop(path=pwd()); Pkg.instantiate()'
julia --project=benchmark/gpu_ace/amd -e 'using Pkg; Pkg.develop(path=pwd()); Pkg.instantiate()'
julia +1.10 --project=benchmark/gpu_ace/legacy -e 'using Pkg; Pkg.instantiate()'
python3 benchmark/gpu_ace/acquire.py --destination /tmp/pfw-public-ace
julia +1.10 --project=benchmark/gpu_ace/legacy benchmark/gpu_ace/export_legacy.jl /tmp/pfw-public-ace/Ti_combined_data_oct24_2023_ace_order_3_degree_20_rcut_6.0_cw_moreV.json /tmp/pfw-public-ace/ti-resolved.toml
julia +1.10 --project=benchmark/gpu_ace/legacy benchmark/gpu_ace/export_legacy.jl /tmp/pfw-public-ace/ti64_combined_data_Apr08_rcut_6.0_order_3_degree_15_custom_reg.json /tmp/pfw-public-ace/tialv-resolved.toml
julia +1.10 --project=benchmark/gpu_ace/legacy benchmark/gpu_ace/oracle_legacy.jl /tmp/pfw-public-ace/Ti_combined_data_oct24_2023_ace_order_3_degree_20_rcut_6.0_cw_moreV.json /tmp/pfw-public-ace/ti-oracle.toml
julia +1.10 --project=benchmark/gpu_ace/legacy benchmark/gpu_ace/oracle_legacy.jl /tmp/pfw-public-ace/ti64_combined_data_Apr08_rcut_6.0_order_3_degree_15_custom_reg.json /tmp/pfw-public-ace/tialv-oracle.toml
julia --project=benchmark/gpu_ace test/gpu_ace_public.jl /tmp/pfw-public-ace
julia --project=benchmark/gpu_ace/amd test/gpu_ace_public.jl /tmp/pfw-public-ace --amdgpu
```

The existing TiAl saved model has a separate pinned producer environment. To
collect its independent fixtures/timings without changing that environment:

```sh
JULIA_LOAD_PATH="$PWD/data/ace_ref:$PWD/benchmark/gpu_ace:$PWD:@stdlib" julia -t4 benchmark/gpu_ace/oracle_tial.jl data/ace_ref/reference_model.json /tmp/pfw-public-ace/tial-oracle.toml
```

The reference model is loaded, never refitted or regenerated. Both oracle scripts
set BLAS to one thread; pass `-t4` to match the recorded four-thread runs. Legacy
oracle timing covers three separate original property calls, whereas TiAl uses
the original combined AtomsCalculators call.

The CPU oracle verification writes a hash-bound `verification-cpu.toml`. Public
benchmark execution refuses an export changed since verification. The producer
never refits, edits source JSON, or substitutes TiAl's radial family. It exports
exact analytic polynomial recurrence/configuration, preserves explicit fitted PI
coefficients, converts complex products algebraically to real harmonics, and checks
radial values/derivatives against the original package. No tabulation approximation
is introduced. The source files lack original independent invariant fitting
coordinates: converted PI rows must not be treated as individually invariant
features. Descriptor selection and invariant fitting explicitly reject them.
Reported paper basis counts differ from the converted density/product/PI counts.
See THIRD_PARTY.md for license, authorship, selected fit and checksum provenance.

```sh
julia --project=benchmark/gpu_ace -t4 benchmark/gpu_ace/run.jl --preset quick --backend cpu --data /tmp/pfw-public-ace --output cpu-quick.toml
julia --project=benchmark/gpu_ace/amd -t4 benchmark/gpu_ace/run.jl --preset quick --backend amd --data /tmp/pfw-public-ace --output amd-quick.toml
julia --project=benchmark/gpu_ace/amd -t4 benchmark/gpu_ace/run.jl --preset full --backend amd --models tial --precision Float32 --strategy compact --samples 10 --seconds 1 --output amd-full.toml
```

Options are `--preset quick|full`, `--backend cpu|amd`, `--models tial,ti,tialv`,
`--data DIR`, `--precision Float32|Float64`, `--strategy auto|staged|baseline|compact`,
`--samples N`, `--seconds S`, `--sizes 4,8,18`, and `--output FILE`.
`--spacing` changes the bcc lattice parameter (default 3.4 Å). Ti–Al–V bulk
labels approximate the alloy atomic fractions; the small oracle cluster includes
all three elements. Sizes are bcc cell repetitions per axis, giving 2g³ atoms. Default quick uses
128 atoms for every available model and 11,664 atoms for TiAl; full extends sizes
to 1,024, 11,664 and 27,648 for each selected model. Large legacy PI cases are
expensive; select bounded sizes explicitly when exploring. Missing public exports
are recorded as omitted, never replaced. CPU-only commands require no GPU and
omit GPU cases explicitly through backend choice. AMD commands require the named
RX 6800 XT and fail if unavailable.

Every case warms and checks numerical predictions before timing. Float64 uses a
common Float64 CPU reference; Float32 also checks the rounded CPU model separately
from application precision differences. `prepared_EFV` times numeric evaluation
with resident model/inputs and excludes extraction and transfer. GPU work is
synchronized, uses borrowed outputs, and includes fresh sorted cells/connectivity,
readout refresh and status readback. The worker warms the public GPU calculator once with the requested strategy;
subsequent timed complete calls retain that mapping (`complete_strategy`).
Complete calls include public extraction,
preparation/upload, execution and owning unitful host results. CPU prepared casts
are outside timing. `evals=1` is explicit; raw samples, median/min/max/spread,
throughput, host allocations, model hashes, commit/dirty state, threads/BLAS,
versions, dimensions and memory payloads are recorded in versioned TOML. Diagnostic
stage times insert synchronization and must not be summed into an undisturbed
throughput estimate. Host allocation counts do not measure device allocations.
Original oracle E/F/V timings are generated once with fixtures in the separate
producer and label the exact three-call boundary.

CPU thread variants must be separate sequential Julia processes. Run no competing
tests or benchmarks when collecting timings. Examples:

```sh
julia --project=benchmark/gpu_ace -t1 benchmark/gpu_ace/audit.jl data/ace_ref/resolved_model_v1.toml /tmp/ace-audit-t1
julia --project=benchmark/gpu_ace -t4 benchmark/gpu_ace/run.jl --backend cpu --models tial --sizes 4 --output cpu-t4.toml
```

The SIMD audit records scalar/SIMD, widths 1/4/8, both precisions and existing
angular/fusion alternatives, with LLVM/native artifacts for representative warmed
signatures. Public analytic custom radials use scalar lanes, documented separately
from imported spline SIMD. The equivariant CPU feature evaluator remains a separate
serial path; no GPU feature or new feature-SIMD claim is made.

Full runs also time a five-step deterministic displacement sequence for the smallest
selected size. AMD frames are uploaded outside timing; timed device-to-device
position copies and evaluations include cached-neighbor reuse and forced rebuilds
with skin up to 0.6 Å, bounded by the single-image cell margin. Cache reset happens in per-sample setup, and each sample replays
the identical sequence. CPU sequence timing includes extraction/fresh topology.
The sequence's extra device buffers and final rebuild count are recorded.

For the complete bounded preset matrix, use the sequential driver:

```sh
julia benchmark/gpu_ace/suite.jl --preset quick --backend both --data /tmp/pfw-public-ace --output quick-results
julia benchmark/gpu_ace/suite.jl --preset full --backend both --data /tmp/pfw-public-ace --samples 5 --seconds 0.3 --output full-results
```

Quick launches one Float64 configuration per available backend. Full adds Float32,
a denser bcc workload, separate CPU thread-count processes, staged/baseline/compact
AMD comparisons, and CPU SIMD audits. It avoids the full cross-product. Override
`--cpu-project`, `--amd-project`, `--threads`, `--models`, `--sizes`, `--audit`,
`--samples` and `--seconds` to bound resources without editing code. The per-case
worker additionally accepts `--parallel true|false` and `--spacing ANGSTROM`.
Both-backend probing records unavailable GPU cases as omitted; explicit AMD jobs
fail on unavailable hardware. Every subprocess/log/report is listed in suite.toml.

Report fields and timing/memory exclusions are specified in [REPORTS.md](REPORTS.md).

Float32 implementation checks compare against the rounded CPU model: force
atol/rtol are 2e-3 / 1e-3, virial 2e-2 / 1e-3, and total-energy error is bounded
by 16eps(Float32) times max(abs(rounded total energy), atom count). The latter
allows reduction-order rounding of large fixed reference offsets. These bounds
are implementation checks for these workloads, separate from descriptive errors
against the original Float64 model and user application tolerances. Reports
record all three rounded-CPU errors.

The measured automatic rule and calibrated ranges are documented in
`docs/gpu_ace.md`. The raw dated results preserve explicit staged/baseline/compact
runs, including slower cases. The CPU SIMD investigation is in
`docs/gpu_ace_cpu_simd_audit.md`.
