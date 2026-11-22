# CPU site chunks

`NeighbourList` owns the experimental chunk executor. `NeighbourInterface`
remains a pure protocol, and ACE supplies adapters to its existing single-site
and batched mathematics. One worker runs synchronously on the caller. Multiple
workers use persistent tasks on Julia's default pool. The caller controls CPU
affinity.

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

A callback handles one chunk. `init_scratch()` runs on the caller once per
participating worker per pass and must provide exclusive mutable scratch.
Use `slot_scratch=true` to receive the 1-based worker slot index and return
caller-owned scratch from a slot cache. The example above uses the ACE
convenience scratch adapter with one worker. For a generic threaded ACE
callback, create an independent `ACEEvaluationState(calc)` per scratch slot
or use the threaded calculator route below. Start with one occupied cell per
chunk in cell order and split
its centers when the actual union of destinations would exceed the configured
capacity. An oversized single site fails during layout construction, before
any callback. Destinations are sorted cell-list indices; contiguous runs merge
into one cell-ordered final array per scatter field. Each scatter vector has at most `scatter_capacity` entries. Geometry,
metadata, one shared final array per scatter field and retained features scale
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
All worker/site/output views are borrowed only during the callback. The executor merges scatter and sums before
reusing them. `get_results` returns only after a complete
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
requires new neighbors/workgroups. Stale environment versions invalidate fields
and reject execution/export. Use `refresh_positions!(workgroups, X)` for new
numeric Å positions in original order with the same atom identities, cell
mapping, neighbor graph and periodic image labels. The caller guarantees those
conditions. Refresh invalidates the latest result and preserves completed
site fields. Owning exports from previous passes remain valid.

## Worker ownership and reduction

`form_numa_workers(neighbors, groups; scatter_capacity=4096)` requests at most
`groups` logical owners. Each owns a contiguous range of whole occupied cells
and their atoms, partitioned approximately by site count. Owners are independent
of hardware and worker count. The plan gives workers contiguous runs of chunks
weighted by site and active-edge counts. Neighbor coordinates and previously
published features are available across all ownership boundaries.

For a threaded pass, workers process successive waves of up to eight chunks per
worker. Each chunk in a wave uses a distinct bounded scatter buffer. Workers
acknowledge completion through buffered channels. Merge jobs route contiguous
scatter runs to disjoint parts of destination ownership ranges, including
subranges within a single group. Each destination has one writer during a merge.
Merge jobs acknowledge completion before the next wave can overwrite the slots.
Force, energy and virial reductions use a fixed order for a given plan. The plan
is cached on the workgroup while the topology, mapping and effective worker count
remain valid.

An exception in a worker callback or merge drains all submitted jobs in that
phase and propagates the worker's exception. A failed pass exposes no result and
keeps previously published site fields. Retry with a fresh pass on the same
workgroup. Beginning a pass, exporting or closing while the workgroup is active
raises `ArgumentError`. `close(workgroups)` wakes and joins the idle task pool and
releases its scatter storage. The next threaded pass starts a new pool. Call
`close` when retiring a workgroup that used threads.

For the AtomsCalculators boundary:

```julia
using PotentialFrameWork.ACE: ACEEvaluationState
import AtomsCalculators as AC
state = ACEEvaluationState(calc; scatter_capacity=4096)
r = AC.calculate((AC.Energy(), AC.Forces(), AC.Virial()), system, calc, ps, state)
@assert r.state === state
```

## Opt-in threading through AtomsCalculators

```julia
props = (AC.Energy(), AC.Forces(), AC.Virial())
state = ACEEvaluationState(calc; scatter_capacity=4096)
try
    single_group = AC.calculate(props, system, calc, ps, state;
                                executor=:threaded, workers=4, groups=1)
    spatial_groups = AC.calculate(props, system, calc, ps, state;
                                  executor=:threaded, workers=4, groups=4)
    # Use only if the moved system has the same graph, mapping and image labels.
    refreshed = AC.calculate(props, moved_system, calc, ps, state;
                             executor=:threaded, workers=4, groups=4, reuse_topology=true)
finally
    close(state)
end
```

`executor=:threaded` defaults to `workers=Threads.nthreads(:default)` and
`groups=1`. One requested worker or a single chunk takes the direct serial path.
The pool starts
lazily and is reused across calls, including fresh layout rebuilds. The number
of workers handling a pass is capped by its chunk count. The evaluation state
owns exclusive slot-local angular and numerical scratch, a preparation cache
populated from the explicit parameters, and the worker pool. Keep parameters and
model data read-only during a pass. Distinct states can share a calculator and
parameters for overlapping calls. An active state rejects overlapping calls
and shutdown. `close(state)` joins its workers, drops the worker scratch and
allows later reuse.

The default call rebuilds geometry. With `reuse_topology=true`, the state must
already have a matching chunk layout. Cell, periodicity, atom count and species
channels must match. The caller also guarantees the same atom identities and
cell mapping, active cutoff edges and periodic image labels. Coordinates are
extracted and evaluated again. Rebuild after any topology or mapping change,
including a change of atom order within one species. This is an explicit
fixed-topology refresh with no skin or automatic validity check.

The high-level `AC.potential_energy`, `AC.forces`, `AC.virial` and combined
methods accept `executor=:threaded` and use the calculator's default state
when no state is supplied. Use `close(calc)` for its worker pool. Explicit
`AC.calculate` returns the supplied state unchanged. `batch_size=1` uses the
single-site kernel and larger existing calculator batch sizes use the batch
callback. Selected centers can contribute forces to any neighbor.

Existing defaults remain available through the explicit `executor=:chunks`
serial route and `executor=:colored` reference route. Legacy
`parallel=true` defaults to coloring for calls without an explicit
`ACEEvaluationState`. Legacy convenience energy-only calls and non-numerical
parameter paths keep the functional evaluator. The explicit-parameter
ChainRules rule can use a threaded primal and continues to build an owning
numerical snapshot for the serial backward. Lux and the numerical derivative
boundary are documented in [site_pullbacks.md](site_pullbacks.md).

For timed breakdowns, pass a fresh `NeighbourList.ChunkTimings()` as `timings`.
The ordinary forward route runs without the stage timers. See the
[workstation report](../benchmark/cpu_chunks/THREADING-REPORT.md) for the
benchmark boundaries and reproduction commands.
