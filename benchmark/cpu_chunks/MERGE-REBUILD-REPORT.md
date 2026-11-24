# CPU ACE force merge and fresh neighbor rebuilds

The starting revision was `3ddd74c8fb90177d31414915ddff285e81c228ae`. The handoff
document was the only untracked file at inspection. The existing reference
fixtures and PR #20 results were preserved. Each raw TOML records hashes of
the loaded source, both harness files, the model and the fixture code. The
merge measurement used the starting source plus the saved merge-only patch.
The rebuild baseline used the unmodified starting source.

All runs used Julia 1.13.4 on a Ryzen 9 9950X with the `powersave` governor and boost enabled,
BLAS=1 and zero interactive Julia threads. CPUs 0–7 are the eight physical
cores in L3/CCD 0. CPUs 8–15 add L3/CCD 1. CPUs 16–31 are the SMT siblings.
Julia tasks can migrate within the process mask. The harness checks the
selected mask against the CPU-to-core mapping and
saves that mapping and the observed governors. The package versions were
AtomsBase 0.5.3, AtomsCalculators 0.2.3, SpheriCart 0.2.4, SIMD 3.8.2,
StaticArrays 1.9.25, Unitful 1.29.0, KernelAbstractions 0.9.48 and
AcceleratedKernels 0.5.13.

Both 3,456-atom BCC fixtures used Float64, SIMD width 4, batch size 1,
scatter capacity 4,096 and one logical owner. Primary rows use ten warmed
public `AC.calculate` calls per variant, cycling the same three frames and
alternating which variant ran first within each pair. These calls export
owning original-order energy, forces and virial. Stage diagnostics ran
separately.

## Force merge

Both policies used the same state, exact directed graph, cell/chunk layout,
waves and routes. The site workers ran in both cases. The coordinator ran
merge routes in the inline policy, with the same per-destination addition
order and serial energy/virial reduction. All three replay frames matched
the colored reference and produced exactly equal E/F/V between policies.
The default remains worker dispatch.

Wall time in ms, median [min–max]. Merge-phase measurements use 11
diagnostic calls per policy. The paired column is the sample-wise public
latency difference, inline minus workers.

| Cores | Fit | Workers public | Inline public | Paired difference | Workers merge phase | Inline merge phase |
| ---: | --- | ---: | ---: | ---: | ---: | ---: |
| 8 | Ti | 131.84 [130.56–137.89] | 131.83 [130.71–134.29] | −0.35 [−5.32–2.85] | 0.170 [0.140–5.772] | 0.176 [0.139–0.226] |
| 8 | Ti–Al–V | 293.48 [292.57–298.95] | 295.49 [291.57–299.86] | −0.08 [−1.89–5.22] | 0.222 [0.188–3.024] | 0.226 [0.182–0.267] |
| 16 | Ti | 104.79 [92.40–115.20] | 85.85 [80.02–95.53] | −15.68 [−25.96–−8.01] | 11.911 [5.723–28.227] | 0.189 [0.125–0.264] |
| 16 | Ti–Al–V | 216.02 [195.35–235.13] | 204.55 [191.30–217.11] | −13.86 [−19.69–8.83] | 10.461 [3.349–14.672] | 0.210 [0.163–0.284] |

At eight cores the public differences were within replay spread. At sixteen
cores the merge-phase ranges were disjoint for both fits. Ti improved in every
public pair. The alloy's public ranges overlap and its paired differences
include a reversal.

The following medians sum selected intervals over waves within a call.
Enqueue starts immediately before `put!`. Entry is the first timestamp in
the worker callback. Receipt wait measures time spent around the coordinating
task's `take!` calls. Each interval may overlap intervals in other columns or
other slots. The raw traces include dispatch start/end, enqueue, entry, worker
finish, receipt start/end and serial-sum start/duration for each wave.

| Cores | Fit | Dispatch (ms) | Sum of per-wave max enqueue→entry (ms) | Sum of worker bodies (ms) | Receipt wait (ms) | Serial E/V, workers / inline (ms) |
| ---: | --- | ---: | ---: | ---: | ---: | ---: |
| 8 | Ti | 0.034 | 0.061 | 0.174 | 0.096 | 0.015 / 0.017 |
| 8 | Ti–Al–V | 0.049 | 0.073 | 0.264 | 0.111 | 0.020 / 0.022 |
| 16 | Ti | 0.041 | 11.745 | 0.217 | 11.250 | 0.010 / 0.014 |
| 16 | Ti–Al–V | 1.771 | 9.510 | 0.240 | 6.610 | 0.012 / 0.023 |

The callback entry and receipt intervals account for much of the observed
wall time at sixteen cores. These timestamps cannot distinguish scheduler
wake-up, task migration and OS preemption inside those intervals.

## Fresh rebuilds

The updated scalar chunk path refills the neighbor graph and spatial layout
through the retained `ACEEvaluationState`. It uses one cell build for directed
traversal and cell order. It stages indices and checked image labels and packs
them by species without building displacement arrays for a query snapshot.
Public `neighbour_list` queries continue to return owning snapshots. The
chunk path skips site coloring. The colored path still prepares coloring when
selected, including when called with the retained state.

Fresh rebuilds are timed on the topology-preserving replay to match the
refresh reference and to isolate the cost of rebuilding each frame. The
harness also displaces an atom across neighbors and verifies both serial
and threaded fresh E/F/V against coloring on that frame and again after
restoring the original geometry. The focused suite independently checks
exact directed indices and image labels against a brute-force snapshot,
including a graph change with the same atom and chunk counts, periodic
crossings, repeated and self images, shrinking and regrowing the graph,
capacity errors, failed rebuild and retry, route and field invalidation,
and owning output lifetimes.

Public E/F/V latency in ms, median [min–max] over ten observations per cell.
The last column compares the median of paired fresh-minus-refresh samples
within each process and shows how much rebuild time changed after accounting
for changes in the refresh reference.

| Cores | Fit | Baseline fresh | Reused storage fresh | Baseline refresh | Updated refresh | Paired rebuild cost, old → new |
| ---: | --- | ---: | ---: | ---: | ---: | ---: |
| 8 | Ti | 151.8 [147.0–174.0] | 122.2 [120.4–126.9] | 115.0 [112.3–125.3] | 114.6 [110.3–135.2] | 37.0 → 9.1 |
| 8 | Ti–Al–V | 321.3 [318.8–336.7] | 290.1 [288.0–297.7] | 282.0 [279.9–290.6] | 279.4 [276.0–285.1] | 39.5 → 11.3 |
| 16 | Ti | 140.1 [133.0–162.6] | 106.3 [96.6–128.0] | 99.5 [92.9–128.6] | 92.1 [86.2–112.7] | 39.1 → 15.1 |
| 16 | Ti–Al–V | 251.0 [241.6–288.7] | 223.7 [196.2–235.0] | 207.5 [183.9–236.6] | 204.3 [189.5–231.1] | 41.4 → 14.7 |

Median allocation in decimal MB and GC wall time in ms per public call. The
raw files retain exact bytes, GC time and wall time for every sample.

| Cores | Fit | Fresh MB, old → new | Fresh GC ms, old → new | Refresh MB, old → new | Refresh GC ms, old → new |
| ---: | --- | ---: | ---: | ---: | ---: |
| 8 | Ti | 107.96 → 12.33 | 9.58 → 0.00 | 2.59 → 2.59 | 0.00 → 0.00 |
| 8 | Ti–Al–V | 110.84 → 15.08 | 9.51 → 0.00 | 5.34 → 5.34 | 0.00 → 0.00 |
| 16 | Ti | 108.24 → 12.61 | 11.40 → 0.00 | 3.00 → 3.00 | 0.00 → 0.00 |
| 16 | Ti–Al–V | 110.72 → 16.22 | 14.10 → 0.00 | 5.35 → 6.62 | 0.00 → 0.00 |

Substage wall time in ms, median of nine probes per cell. Each entry gives
baseline → updated. Baseline environment timings come from an independent
call starting with an empty `FullEnvironment`. The baseline cell-order, layout
and route probes also ran independently and recorded allocation and GC
metrics. Baseline cell order includes a second workspace build. Updated
values come from the instrumented public fresh call. The cell-permutation
probe copies metadata from the shared workspace. Extraction combines
system and workspace extraction. The stages in the table are contained in
neighbor preparation and route setup. They are separate from primary public
timings. Baseline snapshot time includes the original winding/traversal and
materialization together. Updated traversal includes winding and edge staging.
The dash in the coloring row indicates a stage skipped by the fresh chunk path.

| Rebuild substage | Ti, 8 | Ti–Al–V, 8 | Ti, 16 | Ti–Al–V, 16 |
| --- | ---: | ---: | ---: | ---: |
| Element channel mapping | 0.028 → 0.047 | 0.031 → 0.058 | 0.025 → 0.028 | 0.031 → 0.030 |
| Coordinate and workspace extraction | 0.413 → 0.129 | 0.344 → 0.175 | 0.318 → 0.199 | 0.375 → 0.193 |
| Cell sorting, topology and pair schedule | 0.502 → 1.038 | 0.515 → 1.062 | 0.493 → 1.040 | 0.524 → 1.089 |
| Directed traversal and snapshot or staging | 8.156 → 4.424 | 11.311 → 4.528 | 13.067 → 4.494 | 13.308 → 4.421 |
| Species packing and footprints | 9.976 → 2.051 | 11.653 → 2.962 | 12.594 → 2.624 | 13.543 → 3.294 |
| Full site coloring | 1.714 → — | 1.723 → — | 1.608 → — | 1.637 → — |
| Separate cell order or permutation copy | 0.621 → 0.004 | 0.646 → 0.004 | 0.651 → 0.006 | 0.759 → 0.005 |
| Chunk layout | 13.087 → 0.705 | 13.223 → 0.710 | 15.019 → 0.764 | 16.110 → 0.682 |
| Route planning | 2.628 → 2.509 | 2.681 → 2.652 | 2.727 → 2.649 | 2.775 → 2.717 |

Baseline component probes allocated 46.5 MB per chunk layout and 8.30–8.58 MB
per route plan. Route construction was measured separately from workgroup
preparation in both implementations. The updated allocator retains compatible
storage through the workgroup and its worker pool, and routes are built for the
new generation. The refresh rows show process-to-process variation,
especially with both CCDs active. At 16 cores, the alloy allocated about
1.27 MB more on refresh in the updated run. The current probes do not
isolate that additional allocation within the public pass.

## Reproduction

Run from the current repository root. The source archives keep the local
fixtures and the main checkout outside the baseline environment. The patch
restores the exact source used by the merge runs. Use a depot containing the
package versions above for a matched dependency environment.

```bash
root="$PWD"
revision=3ddd74c8fb90177d31414915ddff285e81c228ae
scratch="$(mktemp -d)"
mkdir -p "$scratch"/{baseline-src,merge-src,baseline-env,merge-env,updated-env}
for variant in baseline merge; do
    git archive "$revision" src ext Project.toml |
        tar -x -C "$scratch/$variant-src"
done
git -C "$scratch/merge-src" apply \
    "$root/benchmark/cpu_chunks/results/merge-rebuild/merge-pre-rebuild.patch"
for variant in baseline merge updated; do
    src="$scratch/$variant-src"
    if [ "$variant" = updated ]; then src="$root"; fi
    run_env="$scratch/$variant-env"
    if [ -f "$root/Manifest.toml" ]; then cp "$root/Manifest.toml" "$run_env"; fi
    JULIA_PKG_PRECOMPILE_AUTO=0 julia --startup-file=no --project="$run_env" -e '
        using Pkg, TOML
        Pkg.offline(true)
        src, env = ARGS
        project = TOML.parsefile(joinpath(src, "Project.toml"))
        deps = merge(project["deps"], project["extras"], Dict(project["name"]=>project["uuid"]))
        open(io -> TOML.print(io,Dict("deps"=>deps)),joinpath(env,"Project.toml"),"w")
        Pkg.develop(path=src)
    ' "$src" "$run_env"
done

bench=benchmark/cpu_chunks/merge_rebuild.jl
for cores in 8 16; do
    mask="0-$((cores-1))"
    taskset -c "$mask" julia --startup-file=no --threads="$cores",0 \
        --project="$scratch/baseline-env" "$bench" rebuild "$scratch/baseline-$cores.toml" 10 12
    taskset -c "$mask" julia --startup-file=no --threads="$cores",0 \
        --project="$scratch/merge-env" "$bench" merge "$scratch/merge-$cores.toml" 10 12
    taskset -c "$mask" julia --startup-file=no --threads="$cores",0 \
        --project="$scratch/updated-env" "$bench" rebuild "$scratch/rebuilt-$cores.toml" 10 12
done
julia --startup-file=no --threads=4,0 --project="$scratch/updated-env" -e \
    'using Test; @testset "PotentialFrameWork" begin include("test/runtests.jl"); end'
```

The [raw results and merge-only patch](results/merge-rebuild/) are alongside
the prior benchmark records.
