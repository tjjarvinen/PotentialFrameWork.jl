# Serial CPU site chunks

`NeighbourList` owns the experimental chunk executor. `NeighbourInterface`
remains a pure protocol, and ACE supplies adapters to its existing single-site
and batched mathematics. This milestone runs one workgroup and one reusable
worker slot synchronously. It installs no NUMA affinity or task pipeline.

```julia
using PotentialFrameWork, StaticArrays
using PotentialFrameWork.NeighbourList
using PotentialFrameWork.ACE: create_scratch

# calc is an ACECalculator and system is an AtomsBase system.
# Passing calc supplies its species-channel ordering, precision and SIMD width.
nlist = prepare_neighbors(system, calc)
workgroups = form_numa_workers(nlist, 1; scatter_capacity=4096)
foreach_worker!(workgroups, 1;
    scatter=(forces=SVector{3,Float64},),
    sum_outputs=(energy=Float64, virial=SMatrix{3,3,Float64,9}),
    init_scratch=()->create_scratch(calc),
) do worker
    buffers = get_buffers(worker)
    scratch = get_scratch(worker)
    for site in worker
        evaluate!(buffers, scratch, calc, site)
    end
end
result = get_results(workgroups) # owning numerical Å/eV outputs, original order
```

Use Float32 field types for a Float32 calculator. `prepare_neighbors(system,
cutoff; width=1, precision=Float64)` is the generic alternative, with provider
species channels. Models requiring a different order must supply `channels` and
`species_count`; provider IDs are never implicitly ACE model IDs. The calculator
overload is the naming extension to the proposed cutoff-only preparation.
`create_scratch` refreshes derived readout weights from explicit parameters;
its optional `parameters` and `properties` keywords select authority/requests.
The convenience scratch adapter grows site-local capacity before evaluation.
Do not update parameters during a pass.

A callback handles one chunk, and `init_scratch()` runs once per pass. Centers
are traversed in cell order. Start with one occupied cell per chunk and split
its centers when the actual union of destinations would exceed the configured
capacity. An oversized single site fails during layout construction, before
any callback. Destinations are sorted cell-list indices; contiguous runs merge
into one cell-ordered final array per scatter field. The reusable scatter
vectors occupy at most `scatter_capacity` entries each; there is no full force
array per worker. Geometry, metadata, final arrays and retained features scale
with the system. Scratch for mathematics scales with the maximum site/batch
size, independently of the number of distinct scatter destinations.

`ChunkBuffers` contains exactly `scatter`, `sum_outputs`, and `per_site` named
tuples. Omitted categories are empty. Scatter fields are writable one-based
vectors of the declared value type. Sum fields are writable `Ref`s initialized
with `zero(T)`. Add contributions to both. All active scatter entries and local
sums reset before each callback, and final totals start fresh for each pass.
Per-site fields are writable one-based vectors backed by the cell-ordered result
field; assign every center using `site_ordinal(site)`. Incomplete assignments and
ambiguous names across output categories fail explicitly. Vector-valued sums
remain sums and are not permuted as atom fields.

`central_slot(site)` maps the center to the shared scatter mapping;
`neighbor_slot(site, edge)` maps an active edge, returning zero for inactive
padding. `NeighbourInterface` geometry/species accessors work on these sites.
Original atom identities, species ranges, edge displacements, periodic image
labels and inactive padding are retained in the underlying complete site view.
Repeated images of one atom share a scatter destination while retaining distinct
geometry. ACE uses +g at the center, -g at the neighbor, and virial -d⊗g.
All worker/site/output views are borrowed only during the callback; the executor
merges the buffer before reusing it. `get_results` returns only after a complete
successful pass and exports requested atom fields in original order.

[The runnable two-pass example](../examples/site_chunks.jl) writes `layer1`, then
reads it with `site_inputs=(:layer1,)` using `central_feature(site, Val(:layer1))`
and `neighbor_feature(site, edge, Val(:layer1))`. Inputs are bound once per pass
and are read-only by contract, including mutable model-defined feature values.
Inactive edges have no feature. Multiple images read the same atom feature.
A pass cannot read and write the same field. Completed fields remain in cell
order across passes; `export_field(workgroups, :layer2)` makes an original-order
copy. Only a fully successful pass publishes new fields. Failed passes leave
previously completed fields intact and expose no partial result. Rebuilding
requires new neighbors/workgroups; stale environment versions invalidate fields
and reject execution/export. Do not mutate prepared geometry or topology directly.

For the AtomsCalculators boundary:

```julia
using PotentialFrameWork.ACE: ACEEvaluationState
import AtomsCalculators as AC
state = ACEEvaluationState(calc; scatter_capacity=4096)
r = AC.calculate((AC.Energy(), AC.Forces(), AC.Virial()), system, calc, ps, state)
@assert r.state === state
```

The supplied state owns environment, angular provider, prepared coefficient
cache, numerical scratch and chunk execution storage. Its `model_state` retains
the model evaluation state. Distinct states can share a calculator and read-only
parameters during overlapping calls; a single state has one active invocation.
Every call rebuilds geometry, with no skin reuse. Selected centers retain full
neighborhoods and may write forces to nonselected atoms. `batch_size=1` calls
the site kernel; larger existing calculator batch sizes use the existing batch
callback with bounded outputs consumed before reuse. State scratch capacity is
reused after growth.

High-level E/F/V methods use a calculator-owned default evaluation state, requiring
one active call per calculator. Legacy state values remain accepted and returned
unchanged. The colored reference is available with `executor=:colored` (and
`parallel=true` for its existing threaded executor). The chunk route always runs
serially when explicitly selected, including if `parallel=true` is forwarded.
Legacy convenience calls with `parallel=true` default to the colored threaded
executor; an explicit `ACEEvaluationState` defaults to serial chunks. Legacy energy-only convenience calls
retain their functional path, as do non-Float parameter values and the existing
Lux/ChainRules APIs. An explicit `ACEEvaluationState` uses chunks for energy-only
inference too. No new AD rule is claimed for this mutable executor. Indexing
transpose tests cover per-site center gather and the force-gradient cotangent
`force_bar[center] - force_bar[neighbor]` after permutation to cell order. Full
ACE force-loss pullback wiring through chunks remains a follow-up.

The next milestone needs persistent Base tasks, buffered Channels carrying chunk
assignments/completions, exclusive slots and a merge-before-reuse handshake, plus
failure propagation and pass barriers. Cross-slot reductions and performance
comparison with coloring need measured threaded validation. NUMA placement,
cross-group boundary reduction and optional Dagger execution follow that work.
