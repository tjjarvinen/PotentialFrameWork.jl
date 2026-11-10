# Controlled-affinity repeats

Float64, 1,024 atoms, ten warmed samples per boundary. Four threads use physical CPUs 0–3; sixteen use 0–15. BLAS uses one thread. External modern public calls use the same CPU masks and thread counts.

| Threads | Model | Historical prepared ms | Current reference ms | Optimized ms | Incremental ratio | Total historical ratio | PFW public ms | Modern public ms |
|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 4 | small-tial | 8.542 | 9.176 | 7.752 | 1.18x | 1.10x | 22.398 | 10.080 |
| 4 | medium-ti | 10.392 | 9.844 | 6.794 | 1.45x | 1.53x | 21.279 | 11.946 |
| 4 | large-tial | 281.494 | 259.220 | 55.174 | 4.70x | 5.10x | 74.034 | 90.587 |
| 16 | small-tial | 6.395 | 5.635 | 4.341 | 1.30x | 1.47x | 18.573 | 10.915 |
| 16 | medium-ti | 6.002 | 6.931 | 5.192 | 1.33x | 1.16x | 18.785 | 8.838 |
| 16 | large-tial | 141.253 | 130.794 | 37.513 | 3.49x | 3.77x | 56.056 | 49.659 |

These repeats confirm gains while showing sensitivity to placement and process/compiler state. Historical/source and current-reference ratios are distinct observations; all differences between their executors cannot be attributed to the product primitive. Raw samples/spreads and independent E/F/V errors remain in the TOML files.

Public comparisons include fresh neighbor preparation and output construction on both sides. Prepared PFW timings exclude those costs and are not compared as pure kernel speedups against external public calls.
