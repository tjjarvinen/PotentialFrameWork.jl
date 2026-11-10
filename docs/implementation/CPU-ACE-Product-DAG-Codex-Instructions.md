# CPU ACE product DAG — Codex implementation instructions

Date: 10 November 2026
Repository: https://github.com/tjjarvinen/PotentialFrameWork.jl
Planning reference: `4985499b6eb12421f512c3cb9d27671e5f76ec9b`

## Goal and scope

Implement shared product evaluation and reverse propagation for unrestricted prepared CPU scalar ACE inference. Measure the incremental effect against the completed stage-one prefix/reverse implementation, then stop and produce a report for review. The purpose is to understand performance step by step and collect evidence useful for later MACE work.

This is an implementation task for the executing Codex agent. Proceed through code, validation, measurements and reporting. Do not begin the next optimization stage after writing the report. Do not spawn subagents unless separately authorized. Preserve unrelated work and supported Julia versions.

Keep radial interfaces, SIMD width defaults, neighborhood layout/preparation, scheduling, reductions, GPU execution, functional AD, equivariant feature execution and model/checkpoint formats unchanged. Do not introduce MACE layers, across-site batching, CPU skin reuse, or a PACE exporter. Small integration changes needed to hold a plan and its scratch are in scope. Do not combine the DAG with a contraction or bounds-check optimization campaign.

## 1. Inspect and preserve the baseline

Read repository AGENTS.md, implementation STATUS, and current code before choosing types/signatures. Do not reset the checkout to the planning reference. Record actual starting commit and dirty state; preserve a reproducible source snapshot if dirty. Never overwrite historical results, model caches or pinned manifests.

Relevant files at the planning reference:

- `src/submodules/ACE/projection.jl`: polynomial reference, `unrestricted_products!`, `prepared_products!`, algorithm validation.
- `src/submodules/ACE/cpu.jl`: `SiteScratch`, `PreparedACE`, `CPUExecutor`, `NI.prepare_potential`, `NI.evaluate_site!`, preparation and weight refresh.
- Calculator/interface dispatch forwarding `product_algorithm`.
- `benchmark/cpu_ace/README.md`, `REPORT.md`, drivers, diagnostics, and recorded stage-one reproduction commands.
- Existing product, CPU executor, interface, AD and equivariant tests.

At this reference, `:reference` means polynomial projection and `:optimized` means the independent division-free prefix/reverse algorithm when `policy.orders === nothing`. Preserve those meanings and keep both runnable. Add an explicit `:dag` choice consistently through existing supported entry points. Keep `:auto` on its existing behavior until correctness and measured benefits justify a documented promotion. Restricted-order policies must retain polynomial projection, including when `:dag` is requested, following the existing fallback convention. Report the effective algorithm, not only the requested symbol.

Preserve the stage-one commit as the immutable baseline as well as using same-checkout switches. Shared integration changes could affect all switches.

## 2. Mathematical contract

For unrestricted projection define `x[a] = U[a] + V[a]`. For each original product coordinate j:

    P[j] = product(x[a] for a in products[j])
    Eproducts = sum(weights[j] * P[j] for j)

The DAG must produce P in exactly the existing coordinate order, with weighted density adjoints satisfying the existing dU/dV contract. Both adjoints receive the derivative with respect to x. Current precontracted weights remain authoritative on every execution, including masks and explicit parameter updates. Retain the existing separate energy accumulation for the first implementation so its change does not confound the DAG comparison.

Do not merge fitting coordinates, change coupling matrices, prune fitted zero weights, or specialize the graph to one parameter snapshot. Duplicate algebraic products may share a terminal node, but every original output coordinate must remain mapped and every terminal weight must contribute. Outputs that are also internal nodes receive both their own readout seed and contributions from dependent nodes.

Multiplication and reverse must be division-free. Repeated factors, one or several zero factors, negative values, zero weights, and shared descendants must work. For a multiplication `z = l*r`, reverse uses additive updates:

    adj[l] += adj[z] * value[r]
    adj[r] += adj[z] * value[l]

When l and r refer to the same node both contributions are required. No stale adjoints or output values may survive calls, changed weights, or changes of algorithm. Handle constants/empty products only where supported by the actual model contract; otherwise reject them consistently.

## 3. Build one simple deterministic plan

Construct the graph outside site loops from structural product indices only. Start with a deterministic prefix-sharing trie/DAG over canonical factor sequences (preserve multiplicity). Intern multiplication nodes by their operands, keep leaf density indices distinct from internal node IDs, and store an original-product-to-terminal mapping. Canonical factor sorting is permissible for these scalar commutative products; validate the floating-point reassociation error explicitly.

Use compact topologically ordered arrays for operand indices and terminal mapping. Evaluate forward in topological order and reverse in the opposite order. Plan data is shared read-only by task slots; numerical node values and adjoints belong exclusively to each existing scratch slot. Include a structural compatibility check or reconstruction boundary so trimming, changed topology, or a new resolved model cannot reuse a stale plan. Coefficient changes alone must not rebuild the graph.

Build at model/executor preparation, never per site or per evaluation. Inspect how complete public calls prepare executors: include whatever plan construction they actually incur in public timings. Do not add a global mutable cache to hide this cost. Retain constructor compatibility, scratch growth behavior, and the generic neighbor-interface preparation route.

Keep this first graph strategy bounded. Do not search many heuristic factorization orders or implement a global optimal common-subexpression search. Record how much sharing the simple strategy achieves. If it fails to reduce useful work or is slower, report that result and its evidence; further graph heuristics are a later decision.

## 4. Scratch and integration

Allocate DAG node values and adjoints once per scratch slot; reuse them for every center. Reuse leaf storage where safe, or include any separate leaf arrays in byte counts. Clear the required adjoint region each invocation. Keep polynomial scratch available for fallback, but document the total cost of retaining both paths.

Do not allocate a fresh graph per slot. Preserve slot ownership without relying on `Threads.threadid()`, safe colored scheduling, the single shared force array and owning output lifetimes. Retain zero host allocation for warmed serial prepared evaluation on the existing representative cases.

Update diagnostics to include the real DAG kernel and real scratch arrays. Existing diagnostics hard-coded to the prefix kernel must not be relabelled as DAG measurements. If adding fields to PreparedACE or SiteScratch, update compatibility constructors and every construction/resizing route. Do not serialize transient DAG scratch or change saved model coordinates.

Keep functional AD and equivariant feature paths on their existing algorithms. The DAG reverse here calculates density adjoints for inference forces/virial; it does not establish differentiability of the new mutation path. Record the path used by existing parameter/position derivatives.

## 5. Correctness gates

Compare DAG against both prefix/reverse and polynomial reference, for Float32 and Float64:

- All supported product orders, repeated indices, duplicate products, shared prefixes, an output reused as an intermediate, unused leaves, mixed signs, zero factors and zero weights.
- Nonempty unrestricted selected-species partitions, with meaningful nonzero U and V coverage.
- Repeated calls after changing inputs, readout weights and masks; switching algorithms on the same executor; scratch growth and changed model topology.
- Weighted-adjoint finite differences for a small independent subset, plus direct product values. Use established precision-appropriate tolerances; explain reassociation differences instead of loosening existing contracts.
- Restricted-order policies demonstrably use the general path and match prior results.

Run existing full-model E/F/V comparisons against stage-one PFW and independent cached modern/legacy oracles on perturbed nonsymmetric structures. Retain force/strain finite-difference coverage, masks/trimming, central selection, custom radials, explicit parameter refresh, generic neighbor-interface preparation, units and supported geometries. Exercise serial and parallel results, constructor compatibility and owning output lifetime.

Run the appropriate existing full package/regression suites, including AD, equivariance/parity and CPU-only loading. Test both supported Julia versions used in stage one where available; otherwise report coverage precisely. GPU execution is unchanged; only claim device checks actually run. Preserve established tolerance contracts.

## 6. Controlled benchmark sequence

Reuse existing fixtures and BenchmarkTools machinery. No model refitting or new acquisition is needed. Run expensive benchmarks sequentially; separate compilation from warmed trials. Use at least ten warmed samples for representative comparisons, record every actual count and spread, and flag under-sampled cases. Keep budgets configurable.

First freeze all non-product settings, including precision, width, fusion, angular backend, BLAS threads and affinity. Use the same width for baseline and DAG; do not select a favorable width separately for each.

1. **Serial product diagnostics:** measure forward and reverse separately where feasible, combined product/reverse, and existing readout. Cover small TiAl, medium Ti, large TiAl, published Ti, and published Ti–Al–V, in both precisions. Use representative real site inputs, not only zeros or one synthetic vector. Keep instrumentation out of production timings.
2. **Serial complete evaluation:** compare prefix versus DAG for prepared combined E/F/V at approximately 128 and 1,024 atoms for all fixtures, with representative small and large modern cases at 11,664 atoms. Include both precisions at 1,024 atoms and selected larger cases within the configurable budget.
3. **Public calls:** measure actual complete combined calls, including plan construction/preparation and owning outputs. Separately record one-time plan build latency/allocation and executor preparation. Explain the reuse boundary and, if relevant, how many calls amortize preparation costs.
4. **Thread scaling:** after serial gates, compare representative small/large modern 1,024-atom cases at 1/4/16 physical cores when available, using controlled affinity and BLAS=1. Record fewer available cores honestly. Do not change scheduling in response to these results.
5. **External context:** reuse cached independent oracles for validation. Rerun selected matched modern ACEpotentials public timings if needed for a current speed comparison; label cached measurements and timing-boundary differences. Legacy ACE1 comparisons retain their existing distinct cache/call boundary. PACE stays deferred and nonblocking.

Report regression cases too. Stage microbenchmarks perturb caches and are not additive runtime percentages. Sampling profiles of undisturbed prepared/public workloads should anchor the final bottleneck assessment.

## 7. Structural and memory evidence

For each model record density leaves, original product count/order distribution, unique terminal count, multiplication-node count, intermediate-only nodes, forward/reverse operation estimates for prefix and DAG under clearly stated counting conventions, maximum graph depth, and a sharing ratio. Include clearing, leaf formation and seed accumulation costs rather than equating fewer multiplications with an assured speedup.

Record immutable plan array payload, per-slot numerical DAG scratch, retained fallback scratch, total per-slot scratch, all-slot scratch at measured thread counts, shared force payload and observed host allocation volume. Distinguish these from allocator peaks and hardware cache capacity. Report construction allocations separately from steady execution.

Inspect actual concrete Float32/Float64 DAG callee code for unexpected widening, allocation and indexing costs. Hardware counters/cache measurements are useful only when available and properly collected. If unavailable, say so; do not infer bandwidth saturation from irregular access, static instruction counts, or scratch size alone.

## 8. Report and completion

Save fresh raw machine-readable observations, exact commands and environment/source/model/geometry hashes under a new result directory. Preserve stage-one artifacts. Write a dedicated DAG report linked from the CPU benchmark README and update implementation STATUS with completed checks and limitations.

The report must make these comparisons easy to assess:

- Polynomial reference → stage-one prefix → DAG, with incremental speedup against prefix as the main result.
- Product kernel gains versus prepared E/F/V and public-call gains.
- Model complexity/sharing versus speedup or regression.
- Graph construction and scratch growth versus any evaluation savings.
- Float32/Float64 numerical disagreement and representative external errors.
- Serial versus thread scaling, with effective algorithm and identical settings recorded.

Include a short MACE interpretation: shared scalar products and reverse accumulation may inform later tensor-product execution, but ACE monomial sharing does not establish MACE tensor contraction speedups, training AD speedups, or GPU behavior. State transferable observations and unresolved hypotheses separately.

Completion requires a selectable validated DAG path, preserved prefix/reference/fallback paths, passing relevant regression gates, reproducible timing/memory evidence, and a measured recommendation about the next bottleneck. No fixed speedup is required. Promote `:auto` only for a documented supported regime where repeated complete prepared measurements justify it without hiding public-call regressions. Prefer leaving auto unchanged over introducing a fragile heuristic; explain the decision.

End with the report and a recommendation for review. Do not proceed to radial SIMD, preparation/scheduling, density/force layout, graph factorization experiments, or MACE implementation in the same task.
