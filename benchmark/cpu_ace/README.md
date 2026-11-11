# CPU ACE benchmarks

Measured results and limitations: [stage-one report](REPORT.md) and
[CPU product DAG report](DAG-REPORT.md).

The DAG comparison uses `dag.jl`; its fixed settings, immutable stage-one
baseline, exact sequential commands, raw observations and environment snapshots
are in [DAG reproduction](results/2026-11-10-dag/REPRODUCE.md).

`run.jl` is CPU only. It imports no vendor package and needs no GPU device.
`quick` covers all model/precision combinations at 128 atoms plus small TiAl at
11,664 atoms; `full` defaults to 128/1,024/11,664 atoms for every model. Sizes are
bcc replication counts: `4,8,18`. Model generation is a separate explicit step.
Expensive jobs run sequentially. The driver runs serial first; full defaults to
1/4/16 Julia threads (override to match physical cores). BLAS uses one thread.
Affinity is uncontrolled and reported; use `taskset` externally when appropriate.

## Environments and models

Keep the root Project/Manifest and pinned `data/ace_ref` environment unchanged.
Create a CPU application environment from `benchmark/cpu_ace/Project.toml` and
`Pkg.develop(path=ROOT); Pkg.instantiate()`. Features additionally use the declared
ChainRulesCore/ForwardDiff dependencies. There are no AMDGPU dependencies here.
Modern generation uses an **isolated copy** of `data/ace_ref/{Project,Manifest}.toml`
with BenchmarkTools and StaticArrays added. ACEpotentials must remain 0.10.2;
Julia 1.13.3 was used, while the original pinned reference stays untouched.

```sh
JULIA_DEPOT_PATH=/tmp/pfw-stage1-depot:$HOME/.julia julia --project=/tmp/pfw-stage1-modern -t1 benchmark/cpu_ace/generate_modern.jl benchmark/cpu_ace/results/2026-11-10/models
JULIA_DEPOT_PATH=/tmp/pfw-stage1-depot:$HOME/.julia julia --project=/tmp/pfw-stage1-env benchmark/cpu_ace/suite.jl --project /tmp/pfw-stage1-env --preset quick --output /tmp/pfw-cpu-quick
JULIA_DEPOT_PATH=/tmp/pfw-stage1-depot:$HOME/.julia julia --project=/tmp/pfw-stage1-env benchmark/cpu_ace/suite.jl --project /tmp/pfw-stage1-env --preset full --threads 1,4,16 --seconds 30 --output /tmp/pfw-cpu-full
```

Generation accepts `PFW_CPU_SIZES`, `PFW_CPU_SAMPLES`, `PFW_CPU_SECONDS` (defaults
`4,8,18`, `10`, `3`). Choose 30 seconds for ten samples of slower large cases.
`--samples` and `--seconds` bound each BenchmarkTools trial, not the entire suite.
Actual counts, every observation, spread and under-sampling flags are retained.
The default full grid is intentionally configurable; large-model repeated-frame
samples can take several seconds. Use `--models`, `--sizes`, `--precisions` in the
worker, `--threads` in the driver, and `--widths 1,4,8` for width experiments.
The driver accepts `--baseline-project ENV` for identical historical runs.
See `baselines/README.md` to reconstruct immutable baseline source.

Generation saves exact resolved spline tables, synthetic recipes and independent
modern E/F/V oracles. The supplied small model is never refitted. Medium Ti uses
order 3/degree 10; large TiAl order 4/degree 12. Synthetic weights are deterministic,
nonzero and distinct across central species, with mild index decay. They are
performance fixtures with **no physical accuracy claim**. Nonfinite or excessive
energies/forces fail generation. Modern `ace1_model` here constructs
`ACEpotentials.Models.ace_model`, not the old ACE1 evaluator; concrete backend
and tensor types are retained. The exporter reuses the independently validated
exact spline conversion in `data/ace_ref/export_resolved_model.jl` without running
that script's original model-writing body.

External precision is Float64. PFW Float32 results are checked against both the
rounded PFW reference algorithm and the external Float64 oracle; the latter
measures application sensitivity and is labelled separately. Float32 external
execution is unsupported in this producer, not silently promoted.

## Boundaries and diagnostics

The prepared trial includes combined numerical E/F/V and refreshing product
weights from authoritative parameters. It excludes extraction, graph preparation
and owning output conversion. Public combined calls include actual AtomsBase
extraction, fresh topology/species packing/color scheduling, executor preparation
and unitful owning outputs. Replay measures three perturbed frames, rebuilding
connectivity every frame. There is **no CPU skin reuse**. First warmup timings are
recorded separately; warmed trials use `evals=1` with GC enabled.

`--diagnostics true` (default for serial) separately measures extraction, neighbor
construction, preparation, executor preparation, first-site geometry/radii,
radial/angular evaluation, density, products/reverse, readout, force contraction,
reduction and owning output copies. Raw instrumented preparation observations
retain existing packing/coloring stage timings and total bytes. Stage timings
perturb caches, overlap, and **must not be summed** to reconstruct whole-call time.
Geometry/radii uses an allocating isolated helper; it is not a production kernel.
Scratch array payload and the single shared force-array payload are distinct
from BenchmarkTools host allocation bytes; neither is an allocator peak.

`--codegen DIR` captures optimized LLVM/native text for actual radial callees
(with explicit radial parameters), products/reverse and complete site signatures
in both precisions. `--profile DIR` optionally saves sampled flat prepared-loop
profiles. Code dumps are static evidence, not runtime instruction percentages.
`features.jl OUTPUT [SAMPLES SECONDS]` reports scalar/vector/tensor multi-channel
features and their existing position/parameter pullbacks. These allocating
functional paths preserve general parity and equivariance; scalar inference
optimization does not imply faster AD or features.

## Algorithm selection

`CPUExecutor(...; product_algorithm=:auto)`, `ACECalculator(...; product_algorithm=...)`,
`execute_prepared!(...; product_algorithm=...)`, public calculator combined calls
and `NeighbourInterface.prepare_potential` accept `:reference`, `:optimized`, `:dag`, or
`:auto`. `:reference` always uses the original polynomial algorithm. `:optimized` and `:auto`
use prefix/reverse products; `:dag` uses a deterministic shared product graph.
These unrestricted algorithms apply **only when `policy.orders === nothing`**;
restricted policies always fall back to polynomial projection, even when selected
density happens to be zero. Nonempty unrestricted species partitions are valid.
Weights refresh on every execution. Functional AD, energy-only calculator calls,
and equivariant features keep their existing algorithms and lifetime semantics.
The existing five- and six-argument explicitly typed PreparedACE constructors
remain valid. `ACE.effective_product_algorithm(requested, model.policy.orders)`
reports the effective prepared product path. Auto stays on prefix/reverse.
Graphs depend only on original product structure and density indices, never
weights, and belong to preparation. Prepared model structure is immutable;
changed topology requires a new prepared potential/executor. Each slot owns its
DAG values and adjoints, with polynomial fallback scratch retained. Ordinary
calculator reuse retains the graph; fresh calculator construction builds it.
Functional AD and features continue to use polynomial projection; the mutable
DAG inference reverse supplies density adjoints for forces/virial only.

## Legacy and PACE

Published Ti/Ti–Al–V remain supplementary ACEpotentials **0.6.6 / ACE1 analytic**
cases, never labelled modern. Use `benchmark/gpu_ace/acquire.py`,
`export_legacy.jl`, and the isolated `benchmark/gpu_ace/legacy` environment.
Attribution, CC BY 4.0 and Crown/AWE notices are in
[THIRD_PARTY.md](../gpu_ace/THIRD_PARTY.md). The stage-one legacy worker
`legacy_oracle.jl` emits the same benchmark interchange schema with explicit
three-separate-JuLIP-calls timing labels. It does not export splines.

The worker probes `/opt/lammps-static/bin/lmp -h`, with `--pace` or `PFW_LAMMPS`
overrides, retaining help/version output and executable hash. PACE timing is
explicitly unsupported until equivalent exported spline interpolation, pair
tables and numerical E/F/V conventions are validated. No LAMMPS comparison is
launched by this suite; a probe is not a performance result.

`oracle_modern.jl MODEL_CACHE OUTPUT` refreshes oracle timings from cached model
recipes without rewriting models or earlier results. It defaults to a 60-second
per-trial cap so larger models can reach ten samples. Select the new references
using worker/driver `--oracle-data OUTPUT`. A representative bounded full run is:

```sh
julia --project=/tmp/pfw-stage1-env benchmark/cpu_ace/suite.jl --project /tmp/pfw-stage1-env --preset full --sizes 8 --threads 1,4,16 --seconds 60 --oracle-data benchmark/cpu_ace/results/2026-11-10/oracles-final --output /tmp/pfw-cpu-scaling
```

The driver's `--precisions` controls the precision grid, and `--profile true`
saves serial sampling profiles. The optional `--baseline-project` must point at
an application environment developed against the immutable historical source,
not the modified source with its reference switch. Historical package source
paths and hashes are recorded independently of benchmark harness hashes.

`pinned.jl current|historical MODEL_CACHE ORACLE_CACHE OUTPUT` repeats the
1,024-atom Float64 prepared/public comparisons with caller-controlled CPU
pinning and records `taskset -pc` from inside the process. On this Ryzen host,
CPUs 0–15 are distinct physical cores and CPUs 16–31 their SMT siblings. The
recorded controlled repeats use CPUs 0–3 for four threads and 0–15 for sixteen.
This targets observed between-process variation in the unpinned historical
parallel runs; it adds no production execution policy. Primary serial and raw
unpinned results are retained.

Legacy boundary detail: JuLIP's `neighbourlist(at, rcut)` defaults to
`recompute=false` and reuses `nlist:default` while geometry is unchanged. The
legacy fixed-geometry three-call timings retain this native cache policy; PFW
complete calls rebuild every time. They are supplementary analytic comparisons,
not identical preparation boundaries. The modern public SitePotential method
constructs its default `PairList` per call; its fixed-geometry timing is a fresh
neighbor-list boundary. Original raw labels and oracle observations are retained.
