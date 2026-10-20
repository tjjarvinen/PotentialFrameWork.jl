Run from the repository root after `julia --project -e 'using Pkg; Pkg.instantiate()'`:

```sh
julia --project benchmark/inspect_simd.jl > /tmp/pfw-simd.ll
JULIA_NUM_THREADS=1 julia --project benchmark/run.jl > /tmp/pfw-bench-1.txt
JULIA_NUM_THREADS=4 julia --project benchmark/run.jl > /tmp/pfw-bench-4.txt
PFW_LARGE=1 JULIA_NUM_THREADS=4 julia --project benchmark/run.jl
```

Use a Julia release satisfying Project.toml (validated with 1.10.12). `PFW_REPEATS`
controls the default seven measured samples, after three complete warmups per
case. All public calls rebuild topology and reset output. Fixed-topology timings
use the internal numeric entrypoint and amortise preparation over all samples;
they are not a Verlet reuse claim. Uniform/clustered, single/multiple species,
128/2048 atoms, widths 1/4/8, serial/parallel are covered; optional large case
has 8192 atoms. W=1 uses scalar Float64 arithmetic, W>1 uses SIMD.Vec.

Reports include minimum/median/maximum time, allocated bytes, complete public
calls, system extraction, sorting/SoA packing, candidate adjacency, sparse
schedule construction, species binding, evaluation (fresh geometry plus callback
and immediate accumulation), task partial merge/output permutation, and a
separate geometry-only scan including occupancy counting. That diagnostic is
not an additive decomposition of evaluation. Workspace memory and force/partial
storage, phase occupancy, task work and SIMD lane occupancy are printed. Timings
use Base.@timed; allocation counts include GC-visible Julia allocations, not
process RSS or peak resident memory. No private-force baseline is implemented.

Benchmarks should run without concurrent tests or other benchmarks. CPU affinity
is uncontrolled; results on a shared host have noise. There is no CI speedup
assertion. Inspect `inspect_simd.jl` LLVM for vector arithmetic, and use Julia's
Profile/@allocated on `Execution.execute!` to investigate regressions.

The recorded run is summarised in [SUMMARY.md](SUMMARY.md); raw output is retained
for [one thread](results-1thread.txt) and [four threads](results-4threads.txt).
