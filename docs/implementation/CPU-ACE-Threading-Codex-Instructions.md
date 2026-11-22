# CPU ACE multithreading — Codex implementation instructions

## Goal and execution scope

Implement and validate this update in the local VS Code checkout of [PotentialFrameWork.jl](https://github.com/tjjarvinen/PotentialFrameWork.jl). Run Julia, tests, and benchmarks on the local workstation. The deliverable is a shared-memory threaded forward executor for scalar ACE energy, forces, and virial (E/F/V), with a reusable site-chunk implementation in NeighbourList and measured workstation performance.

Validate on the 16 physical cores of the Ryzen 9 9950X. One shared-memory node is the CPU target. Keep the design suitable for later placement across NUMA domains on a single server with 256+ cores. The GPU target is a single device per simulation with device-specific tiling. This milestone covers CPU forward execution. NUMA binding and memory placement, halo or ghost exchange, MPI, Dagger integration, multi-GPU decomposition, GPU kernel changes, MACE implementation, and threaded pullbacks are outside its scope.

## 1. Inspect the checkout and capture the baseline

1. Read applicable AGENTS.md files and record the current branch, commit, and working-tree status. Preserve local work, manifests, reference models, and prior measurements. PR #19 (merge commit `90d726d08c21aa3e7bbd9b60dc0d17d31369c3f5`) is the planning reference. Implement against the actual local checkout.
2. Read `docs/site_chunks.md`, `docs/site_pullbacks.md`, and the relevant entries of `docs/implementation/STATUS.md`. Inspect the existing NeighbourInterface and NeighbourList chunk and pullback code, the scalar ACE chunk and colored adapters, AtomsCalculators state handling, and the available CPU benchmark fixtures and scripts.
3. Define the benchmark cases and timing boundaries in §5 before changing executor behavior. Capture serial chunk and colored measurements on the starting checkout with those cases. If the harness is developed alongside the update, run the harness against the starting checkout in a separate local worktree or isolated environment as needed. Record the exact source and harness revisions used.

## 2. Ownership and work distribution

Keep these three concepts separate:

| Concept | Responsibility |
| --- | --- |
| Logical spatial group | Owns a disjoint subset of cells and atoms. Placement on CPU threads or memory nodes is a later policy. |
| Chunk | Unit of site evaluation. Contains assigned centers and uses the existing destination mapping, capacity checks, geometry, and batch interface. |
| Worker slot | Owns reusable mutable scratch and bounded scatter storage. A Julia worker task uses it exclusively until its contribution is ready for reduction. |

Extend the existing `prepare_neighbors`, `form_numa_workers`, and `foreach_worker!` path where appropriate. Keep the existing callback and `ChunkBuffers` contracts. Document actual configuration keywords and buffer lifetimes. Keep the execution, scatter routing, and reductions generic in NeighbourList. The ACE adapter should supply model scratch and compute local contributions using the existing site and batched kernels.

Assign compact spatial groups using existing cell data. Start with a simple reproducible partition and inspect its balance using site and active-edge counts. Reuse the ownership and routing plan while the topology and atom mapping stay valid. Rebuild it when either changes. Shared-memory workers read current neighbor positions and read-only feature fields directly across group boundaries. Atom ownership is independent of the number of workers and does not depend on a fixed hardware topology.

Use a bounded set of persistent `Threads.@spawn` workers on the `:default` pool, with buffered Channels or similarly simple Base primitives for repeated pass dispatch and completion. Keep scratch with worker slots and tasks because Julia tasks can migrate. Keep work queues or batched assignments local to a logical group where practical. Avoid a task launch or central dispatcher round trip per small chunk if that becomes a bottleneck. Keep Julia 1.10 compatibility. Avoid blocking communication between iterations of a Base `@threads` loop.

Use a direct serial path for one worker and keep it available for small workloads. Preserve existing selection and defaults on public entry points. Expose a documented opt-in threaded route through the normal scalar CPU ACE calculator path and provide examples for one and several logical groups. Derive any automatic crossover policy from measurements in a later change.

ThreadPinning.jl belongs at the caller or benchmark level. The executor leaves affinity to the caller. OhMyThreads.jl is an optional benchmark comparison if profiling identifies dispatch or load imbalance as significant. To isolate the scheduling effect, compare alternative work distribution with the same spatial groups, chunk batches, scratch, reduction, and outputs. Report any material differences in task or buffer lifetime. Keep the production executor on Base for this milestone.

## 3. Parallel reduction and pass lifetime

- Reuse the existing bounded per-slot scatter footprint and cell-order indices. Route contributions by destination atom ownership, including writes from sites in other groups and repeated periodic images. Reuse contiguous destination runs where possible. The final array is shared, with one writer per destination in each merge phase. Subdivide a large ownership range into disjoint merge ranges so the one-logical-group case can also merge in parallel.
- Define how a worker publishes a completed slot and how all consumers acknowledge it before reuse. A bounded batch of completed slots followed by parallel destination-range reduction is a valid first implementation. Measure its synchronization cost. Size active contribution storage from worker or slot count and scatter capacity. Avoid a full force array per worker or staging every chunk's contributions until the end of the pass.
- Preserve force and virial conventions, per-site assignment semantics, fresh pass totals, and cell-order fields between passes. Reduce energy and virial from partial sums in a defined order for a fixed plan. Use existing numerical tolerances for comparisons across worker counts and partitions. Export owning results in original atom order using the existing permutation.
- State owns its reusable worker resources. Keep the one-active-call guard for the same mutable state. Distinct explicit states using the same read-only calculator and parameters must stay independent. Define and document startup, reuse, and a shutdown path that wakes and joins workers when their pool is retired.
- If a callback, worker, or merge fails, propagate the originating error, drain or cancel the in-flight work safely, and release pass state after all tasks that could write to it are quiescent. Partial outputs or newly written feature fields must never be published. The next call on the state must succeed under the documented recovery behavior. Account for passes with no work and workers waiting on either an input queue or slot acknowledgment. Passes and shutdown must finish without leaving tasks or channels stuck. Borrowed views remain valid only within their callback.
- Keep serial chunk and colored paths usable. Preserve calculator units, output types, state identity, and supported precision and batch modes. Maintain existing AD routing and snapshot lifetime. In particular, verify that differentiable entry points still use their supported path and that an opt-in threaded forward configuration does not advertise a threaded backward rule.

## 4. Local correctness checks

Add focused cases to the existing suites where they expose threading-specific regressions. Run the complete relevant regression suite locally with unchanged numerical tolerances.

- Compare threaded E/F/V against serial chunks and coloring on existing repository fixtures. Exercise one and several logical groups, different worker counts, E/F/V output subsets, supported precisions, and batches of size 1 and greater than 1. Include an actual threaded call through the normal scalar ACE calculator API.
- Force cross-group contributions with several source chunks reaching one destination. Cover selected centers contributing to unselected neighbors, repeated and self periodic images, a nontrivial original-to-cell permutation, inactive padding, forced chunk splitting, and oversized-site capacity handling. Include empty and isolated cases.
- Repeat evaluations with changed coordinates within a valid topology and after a rebuild. Check that old owning outputs remain valid, new passes start clean, and fields tied to stale mappings are invalidated. If the generic threaded path accepts per-site fields, test a second pass that reads center and neighbor fields published by the first.
- Inject a failure while other workers are active. Verify error propagation, absence of partial publication, safe slot reuse, and successful retry. Check overlapping explicit-state calls and the guard on overlapping calls using a single state. Exercise worker shutdown.
- Verify Julia 1.10 compatibility alongside the main local development version. Report exact Julia versions, thread counts, commands, and outcomes. Check that existing pair, GPU, calculator, and AD regression coverage still passes as applicable.

## 5. Workstation benchmarks

Reuse the existing benchmark suite and reference models where possible. Run end-to-end latency benchmarks with internal instrumentation disabled. Reuse case definitions for separate instrumented runs.

| Campaign | Cases and measurements |
| --- | --- |
| Small-system latency | A small molecule and a small periodic cell with compatible ACE fixtures. Run warmed **energy and force** evaluations at **1, 2, and 4 physical cores**. Use a repeatable coordinate sequence that triggers fresh evaluation. Use the supported geometry refresh path and verify that cached results are not being timed. Report absolute wall latency per structure and separate runs for valid neighbor-list reuse and rebuild. If a needed fixture is unavailable, document the actual compatible case used. |
| Large-system strong scaling | Fixed Ti and Ti–Al–V workloads large enough to keep 16 cores busy. For each fixed case, use **1, 2, 4, 8, and 16 physical cores**. Time warmed **E/F/V**, including merging and public output export. Report absolute time, µs/atom/evaluation, speedup relative to the corresponding one-core run, and parallel efficiency. Report preparation and rebuild separately. Label reuse or rebuild according to the path actually exercised. If a baseline lacks an equivalent reuse path, report that difference. |
| Internal breakdown | Time site evaluation, worker idle time and imbalance, time waiting for reusable scatter slots, force merge, and output export. Record pass wall time and the spread across workers. Keep the sum of worker durations separate from wall time and identify stages that overlap. Record allocations and the peak number and capacity of active scatter slots where practical. |
| Baselines | Compare serial chunks, threaded coloring, and the new executor on matching cases. Include ACEpotentials where the same fitted model, precision, outputs, and timing boundary can be matched. Run one and multiple logical groups locally. If testing SMT, label it separately from the physical-core sweep. |

Use one compute thread per physical core in the primary sweep and record CPU affinity, thread counts, Julia and dependency versions, model fixture and precision, system size and neighbor statistics, group and chunk settings, warmup and repetitions, and the neighbor-list policy. Use a controlled BLAS thread count where applicable. Record a representative median and run-to-run spread. Keep initialization or compilation outside warmed measurements. Call a result MD step time only when the measured boundary includes the integrator.

Use published results as context. [Kovács et al., 2023](https://pubs.aip.org/aip/jcp/article/159/4/044118/2904837/Evaluation-of-the-MACE-force-field-architecture) distinguish force latency from full MD throughput. [Fu et al., Symmetrix-XL, §7.3 and Fig. 4](https://arxiv.org/abs/2610.01036) report 1-to-16-thread speedups of 10.9× for complete ASE calls and 12.0× for LAMMPS NVT steps on their 5,000-atom MACE case. Record our workload and timing boundary beside our results. These MACE measurements serve as literature context for the ACE campaign.

## 6. Handoff and completion

Deliver the implementation through the normal local source tree, focused regression coverage, reproducible benchmark commands and raw results, and a short report in the repository's existing documentation layout. Document the public configuration and the group, worker, reduction, and failure contracts. Include the measured effect on small-system latency, comparison with coloring at larger sizes, and any dispatch, imbalance, force merge, or memory bottleneck the timing data exposes.

The milestone is complete when the local correctness checks pass and the workstation report gives before-and-after evidence across the agreed core counts. Roughly **13× speedup at 16 physical cores** on a sufficiently large case is a stretch target. Explain measured shortfalls without changing benchmark boundaries or loosening numerical tolerances. Multiple logical groups on one NUMA node validate ownership and reduction. Physical NUMA locality and scaling to 256+ cores require later measurements on suitable hardware.
