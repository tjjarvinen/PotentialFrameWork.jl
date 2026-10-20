# Recorded benchmark summary — 2026-10-20

These measurements predate the change from compact integer species IDs to
`ChemicalSpecies` in storage/callbacks. They are retained as historical results;
the new species lookup has not been measured with this complete benchmark matrix.

AMD Ryzen 9 9950X (16 cores / 32 logical CPUs), Linux x86_64, Julia 1.10.12,
generic compiler target, Float64, AtomsBase 0.5.3, SIMD 3.8.1, Unitful 1.29.0.
Four-thread process below; serial modes use one worker. Three full warmups and
seven measured samples per case. Affinity uncontrolled. All systems in this table
have mixed species. Public calls include fresh construction; fixed-topology
execution amortises preparation and is an internal diagnostic, with unchanged
geometry, not a public Verlet reuse feature. Times are milliseconds.

| N | Distribution | Width | Mode | Public median [min, max] | Fixed-topology median | Allocated bytes/public call |
|---:|---|---:|---|---:|---:|---:|
| 128 | uniform | 1 | serial | 0.087 [0.082, 0.092] | 0.020 | 138,064 |
| 128 | uniform | 4 | serial | 0.087 [0.085, 0.103] | 0.024 | 138,064 |
| 128 | uniform | 8 | serial | 0.111 [0.107, 0.122] | 0.040 | 138,064 |
| 128 | uniform | 4 | parallel (4) | 0.135 [0.123, 0.147] | 0.077 | 219,280 |
| 128 | uniform | 8 | parallel (4) | 0.184 [0.154, 0.198] | 0.087 | 219,280 |
| 2048 | uniform | 1 | serial | 1.733 [1.720, 1.781] | 0.473 | 2,204,384 |
| 2048 | uniform | 4 | serial | 1.748 [1.742, 1.824] | 0.495 | 2,204,384 |
| 2048 | uniform | 8 | serial | 1.981 [1.965, 2.011] | 0.723 | 2,204,384 |
| 2048 | uniform | 4 | parallel (4) | 1.519 [1.418, 1.535] | 0.219 | 2,269,808 |
| 2048 | uniform | 8 | parallel (4) | 1.662 [1.593, 3.375] | 0.285 | 2,269,808 |
| 2048 | clustered | 1 | serial | 6.027 [5.998, 6.073] | 5.111 | 1,404,720 |
| 2048 | clustered | 4 | serial | 4.555 [4.537, 6.423] | 3.600 | 1,404,720 |
| 2048 | clustered | 8 | serial | 4.489 [4.476, 4.515] | 3.563 | 1,404,720 |
| 2048 | clustered | 4 | parallel (4) | 4.826 [4.681, 5.082] | 3.677 | 1,458,864 |
| 2048 | clustered | 8 | parallel (4) | 4.389 [4.314, 6.318] | 3.457 | 1,458,864 |
| 8192 | uniform | 1 | serial | 7.181 [7.054, 9.221] | 2.006 | 9,183,152 |
| 8192 | uniform | 4 | serial | 7.326 [7.196, 9.244] | 2.101 | 9,183,152 |
| 8192 | uniform | 8 | serial | 8.216 [8.089, 10.106] | 2.985 | 9,183,152 |
| 8192 | uniform | 4 | parallel (4) | 6.291 [6.134, 8.918] | 0.720 | 9,253,088 |
| 8192 | uniform | 8 | parallel (4) | 6.484 [6.140, 9.142] | 0.996 | 9,253,088 |

Raw [one-thread](results-1thread.txt) and [four-thread](results-4threads.txt)
outputs contain all stages, phases, work/lane occupancy, memory and sample ranges.
The one-thread process also measures phased dispatch with only one worker.

Force storage is exactly 3×N×8 bytes: 3,072 / 49,152 / 196,608 bytes for
N=128 / 2048 / 8192. The N=2048 workspace force array is 49,152 bytes in both
one- and four-thread processes. Boundary output is another 3×N result, and cell
energy/virial partials are shared per task, not replicated per worker.

- Uniform mixed N=128, W=4: 36 phases, 7.8% active lanes. tasks_per_phase=[1, 4, 1, 1, 3, 3, 3, 3, 1, 2, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1].
- Uniform mixed N=2048, W=4: 29 phases, 8.5% active lanes. tasks_per_phase=[36, 43, 39, 46, 45, 41, 38, 36, 37, 40, 42, 44, 47, 42, 42, 38, 34, 35, 29, 26, 20, 16, 9, 8, 8, 4, 5, 2, 3].
- Uniform mixed N=8192, W=4: 31 phases, 8.5% active lanes. tasks_per_phase=[184, 178, 177, 190, 182, 178, 170, 170, 179, 181, 173, 179, 175, 164, 163, 152, 141, 129, 102, 90, 66, 47, 37, 27, 20, 18, 13, 8, 3, 1, 1].

Conservative phases and sparse cell populations limit concurrency. Wider vectors
waste lanes on cell tails/inactive candidates; W=8 is not a universal winner.
Clustered cases expose large task-work imbalance; raw per-phase work shows this
directly. Public timings retain serial extraction/topology/schedule costs, so
kernel-only improvements need not translate into complete-call gains. Small
phases can make parallel execution slower. Further tuning should focus on
candidate construction, block packing and phase utilisation, preserving the
validated ownership and write-footprint contracts. No speedup threshold was
imposed, and these measurements do not establish universal performance gains.

Allocated bytes are GC-visible allocations, not peak RSS. No thread-private
force baseline, cache-line alignment experiment, native spill study or NUMA
comparison was performed. See [reproduction](README.md).
