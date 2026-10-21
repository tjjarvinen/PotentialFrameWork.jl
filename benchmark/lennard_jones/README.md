# Argon Lennard-Jones comparison

Compare PotentialFrameWork's scalar/SIMD/phased executor against
[EmpiricalPotentials.jl](https://github.com/JuliaMolSim/EmpiricalPotentials.jl)
0.2.6 through the same `AtomsCalculators.energy_forces_virial` interface.
This directory is a separate benchmark environment. EmpiricalPotentials and
AtomsSystems are only test extras in the main package.

From the repository root, using Julia 1.10.10 or later:

```sh
julia --project=benchmark/lennard_jones -e 'using Pkg; Pkg.develop(path=pwd()); Pkg.instantiate()'
julia --project=benchmark/lennard_jones benchmark/lennard_jones/run_matrix.jl
```

The supplied Manifest records the measured Julia 1.13.2 dependency graph and
uses a relative `../..` path for PotentialFrameWork. A fresh resolution is also
possible from the Project. To retain the relative manifest path during setup, run
`julia --project=. -e 'using Pkg; Pkg.develop(path="../.."); Pkg.instantiate()'`
from this directory instead. On a different Julia minor release, resolve an
environment compatible with that release before running.

Controls (optional):

```sh
PFW_LJ_THREADS=1,2,4,8 PFW_LJ_CELLS=4,6,10,16 PFW_LJ_REPEATS=7 \
  julia --project=benchmark/lennard_jones benchmark/lennard_jones/run_matrix.jl
# Or one process directly:
julia --threads=4 --project=benchmark/lennard_jones benchmark/lennard_jones/worker.jl /tmp/lj-results
```

`run_matrix.jl` starts independent Julia processes sequentially, passing an explicit
`--threads=N` to each. Defaults: 1/2/4/8 threads; FCC cells per side 4/6/10/16
(N=256/864/4000/16384); three warmups and seven timed samples. `PFW_LJ_WARMUP`
controls warmups; at least one warmup/three samples are required. An optional
positional argument chooses the output directory (default `results` here).

`common.jl` constructs deterministic perturbed FCC Argon using
`AtomsSystems.generic_system`: a=5.26 Å, component perturbation up to 0.04 Å,
periodic wrapping, epsilon=0.010323 eV, sigma=3.405 Å, rmin=2^(1/6)*sigma and
cutoff=2.5*sigma. These are illustrative Argon parameters, not a fitted model claim.
The minimum periodic repeat is four, so box lengths exceed twice the cutoff.

**Cutoff matching matters:** EmpiricalPotentials 0.2.6 uses
`U=epsilon*(f(s)-f(sc)-(s-sc)*f'(sc))`, where `f(s)=s^12-2*s^6`, `s=rmin/r` and
`sc=rmin/cutoff`. Its `emin` parameter is `-epsilon`. This tilts in inverse distance,
not the usual distance-based shifted-force LJ. The benchmark-only `ArgonLJ`
callback implements this formula with its analytic derivative; the reference uses
its existing ForwardDiff derivative. Neither reference code nor its assembly is
modified. All atoms must be Ar; real coincident atoms raise an error.

The benchmark aborts unless energies, original-order forces and configurational
virials agree before and after timing (rtol=2e-10, atol=2e-10 in eV/eV per Å).
Tests in `test/lennard_jones.jl`, run by the main package's `Pkg.test()`, also cover
isolated/mixed/full periodic systems, small wrapped grids, scalar and SIMD widths,
cutoff edges and an independent energy derivative. No timing threshold is in CI.

Two regimes are reported separately:

- `rebuild`: complete public calls including system extraction, list construction,
  scheduling, evaluation, reduction and original-order outputs. Framework workspace
  buffers are reused, but topology is rebuilt every call; the reference constructs
  its default PairList every call.
- `fixed`: unchanged-geometry diagnostics. Reference calls receive one prebuilt
  `nlist`. Framework calls use internal `PairPotentials.execute!`/`outputs` with prepared
  coordinates, topology, schedule and bound parameters, then return the same
  static-vector boundary shape. `setup_ms_excluded` records one setup, outside
  timings and amortised over warmups/samples/checks. These are not complete public
  calls or a Verlet reuse claim. PFW setup includes scheduling; EP setup builds
  its directed PairList. There is no shared topology representation between them.

Reference serial runs use `Folds.SequentialEx(), ntasks=1`; threaded runs use
`Folds.ThreadedEx(), ntasks=Threads.nthreads()`. The reference neighbour-list
backend keeps its default CPU configuration: NeighbourLists 0.6 uses sort-based
construction (recorded Julia 1.13 runs); 0.5 uses its legacy serial builder
(available for Julia 1.10 fresh resolutions). Site assembly uses full-size
per-task force partials and directed pair work; this is the external reference's
algorithm. The framework uses unordered pairs and one shared sorted force matrix,
with W=1 scalar and W=4/8 SIMD, serial or phased parallel. Measured differences
therefore concern the entire implementation, including analytic vs AD derivatives,
list ownership/construction, allocations, scheduling and reduction, not SIMD alone.

CSV output records min/median/max milliseconds, allocated Julia bytes, GC time,
absolute energy/force/virial errors, one-time setup, force workspace size, phases,
and framework stage medians. `fixed` reports zero for unused construction stages;
reference stage fields are NaN because its internals are not instrumented. Stage
medians exclude static-vector boundary conversion after output permutation, so
need not sum to the public median. Metadata files contain hardware, compiler
target, thread count, precision and exact dependency versions; stdout records
phase occupancy. Allocations are GC-visible bytes, not peak resident memory.
Compilation, system generation, correctness checks and explicit GC before each
case are outside timings. GC during calls remains measured. CPU affinity is
uncontrolled; do not run concurrent benchmarks/tests when collecting results.


Recorded measurements are in [RESULTS.md](RESULTS.md), with a standalone
[scaling plot](scaling.svg) ([PNG](scaling.png)), CSV files, dependency metadata
and `run.log`. Regenerate the reports after running the matrix:

```sh
python3 benchmark/lennard_jones/summarize.py
# Optional; requires Python matplotlib (not needed to run the Julia benchmark):
MPLCONFIGDIR=/tmp/pfw-matplotlib python3 benchmark/lennard_jones/plot.py
```

The `parallel=false` reference flag selects a sequential **site executor**;
its default CPU neighbour-list builder may still use the process's worker count.
The one-thread process provides the genuinely serial full-call comparison.
`sorted_force_bytes` and `phase_count` are zero in reference CSV rows because
those framework-specific quantities are not instrumented for the reference;
zero does not mean the reference has no force storage.

For a fresh Julia 1.10 resolution, preserve the recorded 1.13 Manifest by moving
it out of this directory, then run setup with Julia 1.10. The Project supports
both NeighbourLists backends; exact measurements here use the supplied 1.13 graph.


Vector-width results are explicit in [RESULTS.md](RESULTS.md): W=1/4/8 serial
comparisons, W=4 versus W=8 ratios with timing ranges, and separate W=8 thread
scaling for rebuild and fixed-topology regimes. Both widths appear in the
[time plot](scaling.svg); the [width-ratio plot](vector-width-scaling.svg) isolates
W=4 time / W=8 time (>1 favours W=8). These reports reuse the recorded matrix,
which already measured and validated W=8 for every size/thread combination;
no new timing samples were substituted.

For Zen 5 / AVX-512 verification, run:

```sh
julia -C native --project=benchmark/lennard_jones benchmark/lennard_jones/inspect_widths.jl
```

The [recorded code inspection](codegen/README.md) includes LLVM and native code
for both the LJ callback and the fused cell traversal at W=1/4/8. On this Ryzen
9950X, Julia 1.13.2 target `znver5`, W=4 callback arithmetic uses YMM registers,
and W=8 uses ZMM registers, e.g. `vmulpd zmm3, zmm0, zmm0`. Eight Float64 lanes
are 512 bits; actual AVX-512 arithmetic is present in the callback and traversal.
The compiler also emits some ZMM arithmetic in W=1/W=4 traversals, so register
presence alone is not a measure of eight-pair processing or throughput.

Using AVX-512 does not guarantee a faster complete call. The measured W=8 path is
slightly slower for the largest system: 6.219 vs 5.749 ms at 8 threads with rebuild,
and 4.941 vs 4.832 ms for fixed topology. At one thread, serial rebuild times are
22.524 vs 22.186 ms (W=8 vs W=4). This data does not isolate the cause; geometry,
masks/tails, endpoint updates and reductions still need profiling. CPU target,
precision, thread scheduling and timings are recorded rather than assuming a
2× gain when widening from four to eight lanes.
