# CPU ACE threading on Ryzen 9 9950X

## Workload and timing boundary

The checkout started on `multithreading_update_v2` at
`90d726d08c21aa3e7bbd9b60dc0d17d31369c3f5` with the supplied threading
instructions untracked. The original source was archived locally before edits.
The baseline uses that archive. Raw records carry the SHA-256 hashes of the
loaded source and benchmark harness, exact versions and process CPU affinity.

The small case with an open boundary is a Ti₂ dimer evaluated with the
elemental Ti fit. The periodic small case is a 16-atom bcc Ti cell. The scaling
cases use the published Ti and Ti–Al–V fits at 3,456 atoms (12×12×12 bcc repetitions).
Both fitted model conversions and their attribution are documented in the
[fixture cache](../cpu_ace/results/2026-11-10/models/README.md). The alloy
geometry uses the existing deterministic approximately Ti–6Al–4V workload.
All primary evaluations use Float64, SIMD width four, site batch size one,
and scatter capacity 4,096. The requested logical group counts are one and four.

Each case cycles through three different position arrays. Before timing, the
harness checks the exact directed graph, image labels and atom permutation
across all frames. It checks E/F/V against fresh coloring and verifies that
requested outputs change along the replay. Warmed observations include
owning, original-order output export. Fresh geometry calls rebuild on every
frame. Chunk reuse calls extract new coordinates and reuse the graph and
routing plan. Prepared-color reuse explicitly supplies the new coordinates
to the existing prepared colored executor and exports owning E/F/V with units.
It uses a lower-level boundary and omits the small `AC.calculate` property and
state wrapping. The starting checkout's serial chunk calculator rebuilds.
On energy-only requests the public colored baseline uses its legacy functional
energy path. The separate prepared-color reuse measurement uses prepared site
mathematics and constructs owning E/F/V containers. All timings are calculator evaluations.

The workstation has 16 physical cores and 32 logical CPUs in one NUMA node.
The affinity mask contains one logical CPU per physical core, with no Julia
interactive threads and one BLAS thread. Tasks can migrate within the mask.
The main campaign uses Julia 1.13.4. The end-to-end runs have stage timing off.
Two warmup calls follow the parity check for each combination. Seven observations per
combination retain wall time, allocations, and GC time in the raw records. Some one-core
measurements overlapped validation on disjoint physical cores. The two- through
16-core timing campaigns ran sequentially after those validations.

## Geometry and reference fit

| Workload | Atoms | Occupied cells | Active directed edges | Packed lanes |
| --- | ---: | ---: | ---: | ---: |
| Ti₂ dimer | 2 | 1 | 2 | 8 |
| Ti periodic | 16 | 1 | 928 | 960 |
| Ti bcc | 3,456 | 216 | 200,448 | 207,360 |
| Ti–Al–V bcc | 3,456 | 216 | 200,448 | 216,760 |

The two large graphs have 58 active edges per atom and 142 coloring phases.
Julia 1.13.4 used LLVM 20.1.8 and target `znver5`. The benchmark environment
used AtomsBase 0.5.3, AtomsCalculators 0.2.3, StaticArrays 1.9.25, Unitful
1.29.0, SIMD 3.8.2, SpheriCart 0.2.4, KernelAbstractions 0.9.48 and
AcceleratedKernels 0.5.13. The host reported the `powersave` frequency
governor during the campaign.

A separate matched comparison uses the supplied fitted Ti–Al model and its
resolved conversion on 16 periodic atoms. Each backend evaluates owning
public E/F/V with fresh neighbors, Float64, Julia 1.13.4 and one physical
core (CPU 15). Across three frames the largest component differences were
(7.28e−12) eV in energy, (2.04e−13) eV/Å in force and
(1.41e−12) eV in virial.

| Evaluator | Median (ms) | Min–max (ms) | Allocated bytes, first sample |
| --- | ---: | ---: | ---: |
| ACEpotentials 0.10.2 | 0.342 | 0.329–3.867 | 746,968 |
| PotentialFrameWork, one worker | 0.448 | 0.437–4.126 | 328,808 |

Each matched row has nine warmed samples. The reference comparison uses its
own supplied fit because the original upstream inputs for the published Ti
and Ti–Al–V conversions were unavailable in the local checkout. Their prior
benchmark provenance is in the [CPU ACE reports](../cpu_ace/README.md).

## Small-system latency

Median wall time in milliseconds per structure. Both small graphs have a single
chunk, so one and four requested groups use the same direct path. The raw files
retain min/max and all samples for both requests and group settings.

| Case | Property | Cores | Serial fresh | Colored fresh | Threaded fresh | Threaded reuse |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Ti₂ dimer | energy | 1 | 0.355 | 1.113 | 0.386 | 0.346 |
| Ti₂ dimer | energy | 2 | 0.374 | 1.108 | 0.356 | 0.349 |
| Ti₂ dimer | energy | 4 | 0.348 | 1.076 | 0.442 | 0.448 |
| Ti₂ dimer | forces | 1 | 0.344 | 0.872 | 0.353 | 0.336 |
| Ti₂ dimer | forces | 2 | 0.365 | 0.895 | 0.361 | 0.362 |
| Ti₂ dimer | forces | 4 | 0.352 | 0.878 | 0.446 | 0.403 |
| Ti 16-atom periodic | energy | 1 | 2.931 | 9.979 | 2.980 | 2.908 |
| Ti 16-atom periodic | energy | 2 | 2.928 | 10.486 | 2.946 | 2.881 |
| Ti 16-atom periodic | energy | 4 | 2.898 | 10.445 | 2.912 | 2.839 |
| Ti 16-atom periodic | forces | 1 | 3.277 | 3.926 | 3.269 | 3.206 |
| Ti 16-atom periodic | forces | 2 | 3.280 | 3.925 | 3.280 | 3.215 |
| Ti 16-atom periodic | forces | 4 | 3.251 | 3.932 | 3.259 | 3.182 |

On one core the new fresh path was within 3% of serial force latency on both
small fixtures. Reuse saved 17 µs for the dimer and 63 µs for the periodic cell
relative to the new fresh path. A single chunk offers no site parallelism.
Raw spreads show more variability in the four-core dimer timings.

## Strong scaling

All calls return owning E/F/V. Each latency cell gives median milliseconds
with the minimum and maximum observation in brackets.

### Ti bcc

| Cores | Baseline serial fresh | Baseline colored fresh | Prepared color reuse | Threaded fresh, 1 group | Threaded reuse, 1 group | Threaded fresh, 4 groups | Threaded reuse, 4 groups |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 720.943 [716.879–724.221] | 708.355 [704.962–746.729] | 710.245 [709.125–714.194] | 752.867 [741.132–753.437] | 719.832 [717.932–735.340] | 751.019 [747.343–753.583] | 718.223 [716.999–719.100] |
| 2 | 716.353 [714.543–719.733] | 373.055 [369.633–375.831] | 363.165 [362.481–370.962] | 403.711 [397.161–404.312] | 368.812 [367.986–369.737] | 404.854 [401.427–405.410] | 369.736 [368.468–381.163] |
| 4 | 720.916 [713.393–738.246] | 209.676 [204.773–211.787] | 194.027 [193.720–194.506] | 240.515 [233.308–250.527] | 206.492 [203.899–213.379] | 240.847 [236.325–242.541] | 205.517 [205.149–206.571] |
| 8 | 802.846 [798.004–806.871] | 144.694 [139.474–147.014] | 113.861 [112.913–115.215] | 153.580 [144.137–156.738] | 116.069 [114.551–124.637] | 151.525 [144.498–155.711] | 117.390 [115.033–125.716] |
| 16 | 807.577 [800.185–809.831] | 124.788 [112.870–171.850] | 121.968 [95.675–155.111] | 144.222 [127.053–152.332] | 95.342 [89.778–113.121] | 137.747 [130.292–146.583] | 99.494 [92.475–110.100] |

Speedup and efficiency use the matching one-core configuration. Times per
atom are in µs for the one-group fresh and reuse boundaries.

| Cores | Color fresh | Threaded fresh, 1 group | Threaded reuse, 1 group | Threaded fresh, 4 groups | Threaded reuse, 4 groups | 1 group µs/atom fresh / reuse |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 1.00× / 100% | 1.00× / 100% | 1.00× / 100% | 1.00× / 100% | 1.00× / 100% | 217.8 / 208.3 |
| 2 | 1.90× / 95% | 1.86× / 93% | 1.95× / 98% | 1.86× / 93% | 1.94× / 97% | 116.8 / 106.7 |
| 4 | 3.38× / 84% | 3.13× / 78% | 3.49× / 87% | 3.12× / 78% | 3.49× / 87% | 69.6 / 59.7 |
| 8 | 4.90× / 61% | 4.90× / 61% | 6.20× / 78% | 4.96× / 62% | 6.12× / 76% | 44.4 / 33.6 |
| 16 | 5.68× / 35% | 5.22× / 33% | 7.55× / 47% | 5.45× / 34% | 7.22× / 45% | 41.7 / 27.6 |

### Ti–Al–V bcc

| Cores | Baseline serial fresh | Baseline colored fresh | Prepared color reuse | Threaded fresh, 1 group | Threaded reuse, 1 group | Threaded fresh, 4 groups | Threaded reuse, 4 groups |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 1876.035 [1867.336–1887.805] | 1931.266 [1922.954–1937.501] | 1742.838 [1735.077–1776.850] | 1745.592 [1740.705–1763.373] | 1716.518 [1713.955–1723.266] | 1739.953 [1735.133–1741.045] | 1704.748 [1703.675–1707.787] |
| 2 | 1777.901 [1773.376–1779.035] | 907.637 [902.904–910.375] | 922.464 [920.815–923.246] | 983.433 [979.488–989.127] | 946.378 [941.444–948.908] | 968.348 [964.740–977.464] | 929.897 [926.046–932.386] |
| 4 | 1746.157 [1733.984–1756.308] | 484.491 [478.842–487.044] | 482.223 [480.866–489.240] | 542.474 [536.456–557.038] | 505.134 [500.892–508.189] | 544.493 [535.479–546.676] | 508.118 [503.548–532.516] |
| 8 | 1725.174 [1717.828–1733.013] | 293.146 [291.631–304.620] | 279.866 [275.595–283.676] | 323.264 [314.065–326.204] | 286.712 [284.221–301.185] | 323.947 [314.768–346.679] | 287.809 [285.354–297.289] |
| 16 | 1739.055 [1726.481–1753.419] | 292.425 [278.938–320.245] | 291.894 [259.750–325.339] | 246.420 [228.676–254.305] | 200.290 [195.172–215.220] | 255.335 [234.977–268.124] | 212.575 [206.975–223.384] |

Speedup and efficiency use the matching one-core configuration. Times per
atom are in µs for the one-group fresh and reuse boundaries.

| Cores | Color fresh | Threaded fresh, 1 group | Threaded reuse, 1 group | Threaded fresh, 4 groups | Threaded reuse, 4 groups | 1 group µs/atom fresh / reuse |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 1.00× / 100% | 1.00× / 100% | 1.00× / 100% | 1.00× / 100% | 1.00× / 100% | 505.1 / 496.7 |
| 2 | 2.13× / 106% | 1.77× / 89% | 1.81× / 91% | 1.80× / 90% | 1.83× / 92% | 284.6 / 273.8 |
| 4 | 3.99× / 100% | 3.22× / 80% | 3.40× / 85% | 3.20× / 80% | 3.36× / 84% | 157.0 / 146.2 |
| 8 | 6.59× / 82% | 5.40× / 67% | 5.99× / 75% | 5.37× / 67% | 5.92× / 74% | 93.5 / 83.0 |
| 16 | 6.60× / 41% | 7.08× / 44% | 8.57× / 54% | 6.81× / 43% | 8.02× / 50% | 71.3 / 58.0 |

At 16 cores with one logical group and fixed-topology reuse, the Ti case ran in
95.3 ms and the Ti–Al–V case in 200.3 ms. These are 7.55× and 8.57× over the
corresponding one-core runs. Compared with prepared coloring at 16 cores,
latency fell by 22% for Ti and 31% for Ti–Al–V. With a rebuild on each call,
the new Ti executor took 144.2 ms versus 124.8 ms for public coloring, and the
alloy took 246.4 ms versus 292.4 ms. The four-group layouts produce similar
times at lower core counts and add 4–6% reuse latency at 16 cores. Group
ownership and task placement were exercised on one physical NUMA node.

At 16 cores, independent graph rebuild, cell ordering, and unthreaded chunk
layout took respectively 20.1, 0.55 and 12.6 ms for Ti, and 16.9, 0.50 and
9.5 ms for Ti–Al–V. Those probes run separately. The threaded fresh path also
builds routes, refreshes prepared parameters, computes sites and exports.
In the first uninstrumented samples for the 16-core, one-group configuration,
the Ti rebuild and reuse allocated about 108.3 MB and 3.0 MB, and the alloy
rebuild and reuse allocated 110.7 MB and 5.3 MB. These are total host
allocations over a call.

## Worker and merge timings

These are separate warmed 16-core reuse observations with timers enabled.
Worker entries are accumulated elapsed intervals for individual slots.
They overlap across slots. The site phase includes dispatch and waiting for
all workers in each wave. Slot wait starts when a worker publishes its
contribution and ends after all merge acknowledgments. It overlaps the
site tail and merge phase.

| Fit | Groups | Public wall (ms) | Site phase (ms) | Site time across slots, min–max (ms) | Merge phase (ms) | Longest slot wait (ms) | Export (ms) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Ti | 1 | 100.85 | 89.35 | 60.07–84.31 | 10.84 | 40.10 | 0.042 |
| Ti | 4 | 94.85 | 88.53 | 58.23–82.94 | 5.75 | 36.03 | 0.035 |
| Ti–Al–V | 1 | 191.83 | 174.00 | 137.24–167.18 | 11.33 | 48.05 | 0.040 |
| Ti–Al–V | 4 | 217.21 | 200.99 | 136.93–182.98 | 14.34 | 78.34 | 0.046 |

For one group, the 16-core site-phase time outside each worker's own site
interval ranged from 5.04–29.28 ms for Ti and 6.83–36.76 ms for the alloy.
The sum of merge-worker intervals was 0.21 and 0.29 ms respectively, while
the merge phases took about 11 ms wall time in these samples. At four
cores the corresponding merge phases took 0.16 and 0.35 ms. These gaps
reflect coordination and scheduling time in the measured phases.

The 13× stretch target was missed. From the measured one-core reuse latencies,
13× would allow about 55 ms for Ti and 132 ms for the alloy. The 16-core
site phases alone exceeded those budgets. The worker spreads and slot waits
show idle time at the wave barriers, with additional coordination visible
during the merge. Hardware-level site-kernel costs remain unmeasured.

The four ownership groups contain 866, 871, 863 and 856 sites, with 50,228,
50,518, 50,054 and 49,648 active directed edges in both fits. At 16 cores the
216 chunks use two waves with 128 peak active buffers, each with capacity
243 force triplets. The force scatter payload is 746,496 bytes across those
buffers, plus one 82,944-byte cell-order final force array. At four cores
there are 32 active buffers with the same capacity and 186,624 bytes of
force scatter payload. Metadata, numerical scratch, and owning exports are
additional. The raw diagnostics include worker intervals, phase walls,
public allocations, active slot counts, and both geometry policies.

## Reproduction

Run from the package root. The main sweep uses the local resolved fixtures
already tracked in the repository. Julia environments and the archived baseline
can live outside the checkout.

```sh
root=$PWD
baseline=$(mktemp -d)
git archive 90d726d08c21aa3e7bbd9b60dc0d17d31369c3f5 src ext Project.toml |
    tar -x -C "$baseline"
export JULIA_DEPOT_PATH="/tmp/pfw-thread-depot:$HOME/.julia"
export JULIA_PKG_PRECOMPILE_AUTO=0
```

Use separate application environments. For both Julia versions, use a copy of
the root `[deps]` and `[extras]` merged into `[deps]`, with the corresponding
`[compat]` entries. Use `Pkg.develop(path=root)` for the working source and
`Pkg.develop(path=baseline)` for the archive in a copy of the Julia 1.13.4
application environment. The resolved package versions in the raw files identify
the actual benchmark set.

```sh
# JULIA points to the Julia 1.13.4 executable.
# BASE_ENV points to the baseline environment. RUN_ENV points to the updated one.
for cores in 1 2 4 8 16; do
    cpus="0-$((cores - 1))"
    taskset -c "$cpus" "$JULIA" --startup-file=no --threads="$cores,0" --project="$BASE_ENV" \
        benchmark/cpu_chunks/threading.jl baseline "benchmark/cpu_chunks/results/baseline-$cores.toml"
    taskset -c "$cpus" "$JULIA" --startup-file=no --threads="$cores,0" --project="$RUN_ENV" \
        benchmark/cpu_chunks/threading.jl updated "benchmark/cpu_chunks/results/updated-$cores.toml"
done
for cores in 1 4 16; do
    taskset -c "0-$((cores - 1))" "$JULIA" --startup-file=no --threads="$cores,0" --project="$RUN_ENV" \
        benchmark/cpu_chunks/threading.jl diagnostics "benchmark/cpu_chunks/results/diagnostics-$cores.toml"
done
```

For the matched fit, set `REF_ENV` to an isolated copy of `data/ace_ref`.
Each file records the numeric outputs for all three frames. Compare those
outputs before comparing the latency.

```sh
taskset -c 15 "$JULIA" --startup-file=no --threads=1,0 --project="$REF_ENV" \
    benchmark/cpu_chunks/threading_reference.jl reference benchmark/cpu_chunks/results/matched-reference.toml
taskset -c 15 "$JULIA" --startup-file=no --threads=1,0 --project="$RUN_ENV" \
    benchmark/cpu_chunks/threading_reference.jl pfw benchmark/cpu_chunks/results/matched-pfw.toml
```

## Literature context

Kovács et al. distinguish small-structure force latency, which excludes
neighbor construction and other MD work, from trajectory performance.
[Evaluation of the MACE force field architecture](https://arxiv.org/abs/2305.14247).

Fu et al. report 10.9× for complete ASE calls and 12.0× for complete LAMMPS NVT
steps across 1–16 physical cores on a 5,000-atom MACE-OMAT-0 case. Their CPU is a
Ryzen 9 9950X3D2 and the optimized paths use retained neighbors. The measurements
in this report use scalar ACE, a Ryzen 9 9950X, 3,456 atoms for scaling, and the calculator
boundary described above.
[Symmetrix-XL §7.3 and Fig. 4](https://arxiv.org/html/2610.01036v1#S7.SS3).
