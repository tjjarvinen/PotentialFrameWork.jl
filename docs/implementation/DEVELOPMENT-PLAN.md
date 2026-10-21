# Development order and acceptance gates

Complete the initial implementation through milestone 5. Milestone 6 is follow-up
planning, not permission to expand this implementation into ACE/MACE.

## 0 — Inspect and settle conventions

Inspect the repository and the current NeighbourLists/AtomsBase/SIMD APIs. Preserve
the existing Project.toml identity and Julia compatibility. Record dependencies,
list ownership, units, displacement/force/virial conventions and supported cells.
Choose whether to adapt NeighbourLists or implement the required cell storage.

Gate: a short STATUS.md decision record, consistent public API sketch and a
package that loads. Do not claim repository APIs were verified if inaccessible.

## 1 — Serial scalar complete force path

Implement extraction, spatial permutation, cells, candidate interactions and
original-order outputs. Add a callable smooth finite-cutoff pair potential and
mixed-species parameters. Implement energy, force and virial plus energy-only use.
Build an independent O(N²) oracle in tests using original atom order.

Gate: numerical reference, periodic, species and permutation tests pass.

## 2 — Explicit SIMD blocks

Implement 1 × W blocks and immediate result accumulation. Include safe masking,
partially occupied blocks and selected widths/precisions. Avoid a global array of
pair-force results. Add a runnable public example showing a user-defined callback.

Gate: scalar/SIMD agreement and tail-mask correctness. Inspect representative
generated code for vector arithmetic; distinguish it from the scalar baseline.

## 3 — One-array parallel execution

Implement cell task ownership, conservative write footprints, sparse conflict
construction and a validated phase schedule. Keep each cell task internally serial.
Use one shared sorted force array, with only small task/workspace partials for
energy/virial. Ensure repeated calls clear/reset all active outputs correctly.

Gate: footprint-disjointness tests and repeated multithreaded numerical agreement,
including wrapped small grids and irregular cell populations. Check actual memory
allocation/scaling; no hidden N-atoms-by-N-threads force array.

## 4 — Reuse and performance measurements

Reuse workspaces and separate geometry refresh from immutable metadata. A skin is
optional; rebuilding each evaluation is acceptable if explicitly documented and
included in timings. Benchmark scalar, SIMD serial and phased parallel execution.
If practical include a thread-private baseline with selective touched-block merge.

Gate: report build, refresh/packing, evaluation/accumulation, schedule and final
permutation costs, plus complete evaluation time. Report allocations, memory,
phase occupancy and variation across runs. No required speedup threshold: explain
regressions and optimise measured bottlenecks within scope.

## 5 — Package usability and handoff

Document AtomsBase input, units, output conventions, callback contract, limits,
threaded usage and benchmark reproduction. Follow repository CI conventions and
run package tests on supported Julia. Add AtomsCalculators integration only if
its current interface can be implemented cleanly; it is desirable but not a reason
to expand scope or invent an API. Otherwise describe the adapter as follow-up.

Gate: an external user can implement a simple potential without writing neighbour
traversal or force reduction. STATUS.md lists completed gates, evidence and limits.

## 6 — Later work

Measure whether fewer shared force domains improve phase utilisation/NUMA traffic.
Then implement species-grouped site views and a minimal basis projection workload,
followed by ACE kernels and model derivatives. MACE/training follow as separate
stages with their own accumulation design.
