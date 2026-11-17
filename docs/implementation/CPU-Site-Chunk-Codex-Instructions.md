# Codex instructions for the serial CPU site-chunk interface

## Goal and scope

Implement the first milestone in [PotentialFrameWork.jl](https://github.com/tjjarvinen/PotentialFrameWork.jl) in the local VS Code checkout. Run Julia tests locally. Start by inspecting the current checkout, applicable AGENTS.md files, and the existing NeighbourInterface, NeighbourList and ACE CPU/AtomsCalculators paths. Reuse existing neighbor geometry, site views, site and batched evaluators, scratch hooks, and force/virial conventions.

Deliver a generic serial chunk executor, integrate CPU ACE E/F/V, and demonstrate a two-pass site feature calculation. Use one workgroup and one worker. Keep the existing executor available as the parity reference. The subsequent milestone is a single-workgroup pipeline using persistent Base tasks and buffered Channels. NUMA placement and cross-group boundary reduction follow measured single-workgroup results. Dagger is a later optional executor. Keep this work in PFW and preserve existing pair and GPU behavior.

## Caller contract

Use this call shape in PFW. Adapt preparation names to the existing API where appropriate and document any naming differences.

```julia
nlist = prepare_neighbors(system, cutoff)
workgroups = form_numa_workers(nlist, 1)

foreach_worker!(
    workgroups, 1;
    scatter = (forces = SVector{3,Float64},),
    sum_outputs = (energy = Float64, virial = SMatrix{3,3,Float64,9}),
    init_scratch = () -> create_scratch(calc),
) do worker
    buffers = get_buffers(worker)
    scratch = get_scratch(worker)
    for site in worker
        evaluate!(buffers, scratch, calc, site)
    end
end

result = get_results(workgroups)
```

Implement traversal and reduction once in the executor. The potential provides scratch creation and local evaluation. Reuse the existing ACE evaluator through an adapter that routes its outputs. Keep both single-site and existing batched kernels usable. A callback handles one chunk. A reusable worker slot supplies exclusive scratch and a buffer across chunks. Borrowed views are valid during the callback only. The executor merges the completed buffer before reusing it. Return from a pass only after its outputs are complete.

## Chunk and buffer contract

```julia
struct ChunkBuffers{Scatter<:NamedTuple,Sums<:NamedTuple,Sites<:NamedTuple}
    scatter::Scatter
    sum_outputs::Sums
    per_site::Sites
end
```

Each tuple contains exactly its requested field names. An omitted category is an empty named tuple. Reject ambiguous result names within a pass.

| Category | Callback storage for declared value type T | Write rule |
| --- | --- | --- |
| scatter | Writable one-based AbstractVector{T}, indexed by chunk-local destination slot | Add by destination atom. All scatter fields share the mapping. |
| sum_outputs | Writable Ref{T} | Add contributions to a chunk total. Merge totals with +. |
| per_site | Writable one-based AbstractVector{T}, indexed by center ordinal within the chunk | Assign one value per center. Back it with the cell-ordered result field. |

Clear active scatter slots and sum accumulators to zero(T) per chunk. Start pass totals fresh. Per-site storage may be uninitialized before the callback. Export per-atom outputs in original order through the existing atom permutation. Keep intermediates in cell order between passes.

- Start with one occupied cell per chunk. Split by centers when the actual union of destinations would exceed the configured scatter capacity. Report an explicit capacity error if a single site exceeds the cap. Validate chunk layout before evaluation where practical.
- A site exposes its ordinal within the chunk, its center scatter slot, and the slot for each active neighbor edge. Preserve edge-specific displacements and image shifts, species ordering, and inactive padding. Distinct periodic edges that reach the same atom map to one destination slot.
- Order destinations by cell-list atom index and merge contiguous runs where available. Handle periodic wrap and sparse runs from actual indices. Use one cell-ordered final array per scatter field. Bound reusable temporary scatter storage by slot count and configured capacity. Requested per-atom outputs, geometry, and neighbor metadata scale with system size.
- Preserve the existing force convention. For displacement r_neighbor - r_center and gradient g, add g at the center and subtract g at the neighbor. Reuse the existing virial convention, including for periodic edges.

## Per-site fields across passes

Store completed named per-site fields in cell order in the execution state. A later call declares site_inputs = (:layer1,) and may write per_site = (layer2 = FeatureType,). Provide center and neighbor accessors such as central_feature(site, Val(:layer1)) and neighbor_feature(site, edge, Val(:layer1)). Multiple images of one atom read the same feature.

Bind inputs once per pass. Treat them as read-only and reject reading and writing the same named field in a pass. Reject missing or incomplete inputs. Publish newly written fields only after successful completion. Invalidate fields when rebuilding changes their atom mapping. Keep feature representation, channel layout and parity model-defined. Demonstrate explicit export of a completed field in original atom order.

## AtomsCalculators and pullbacks

Route explicit calculate(..., calc, ps, st) through the supplied ACE evaluation state. Put mutable environment, executor workspace, and required model evaluation state there. Return the updated state using the existing AtomsCalculators convention. Two calls can share the calculator and read-only parameters when they use distinct states. High-level convenience methods can use a calculator-owned default state with one active call per calculator. Preserve existing supported entry points and differentiation paths.

Validate the indexing transposes with a small adjoint identity check. A per-site cotangent gathers by center. For force scatter, bring output cotangents into cell order and use force_bar[center] - force_bar[neighbor] as the cotangent of each displacement gradient. Keep full ACE force-loss pullback wiring for a follow-up. The chunk contract should also allow position and feature cotangents to be scattered with the same destination mapping.

## Acceptance checks

1. Compare serial chunk ACE E/F/V with the existing executor on repository fixtures using existing precision-appropriate tolerances. Check original-order force output with a nontrivial atom permutation.
2. Cover periodic images sharing one destination, a forced split within a cell, and an oversized single site. Verify the observed scatter capacity stays within the configured limit.
3. Exercise two feature passes. The second reads both center and neighbor features from the first. Verify cell-order access and original-order export against a direct reference. Include a pass with per-site outputs only.
4. Check output reset and field lifetime across passes and state reuse. Compare two overlapping explicit-state calls on the same calculator with distinct states and geometries against serial references. Each invocation uses the serial chunk executor internally.
5. Run the relevant existing CPU tests locally. Report changed entry points, tests and results, any unresolved issues, and what is needed for the Channel executor. Performance comparison against site coloring belongs to the threaded milestone.
