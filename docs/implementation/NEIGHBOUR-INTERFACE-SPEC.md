# CPU neighborhood evaluation interface specification

## 1. Goal and boundaries

Make PotentialFrameWork a practical testing ground for an optional neighborhood
protocol that might eventually live in AtomsCalculators. Implement now; do not
stabilize an ecosystem API or submit a proposal.

Separate:
- Protocol: callback/query declarations and mathematical/data contracts.
- Neighbor execution: topology, geometry, views, scheduling, buffers, reduction.
- Potential mathematics: pair/ACE model structure, explicit parameters, model
  preparation and local scratch.

Put the experimental protocol in a local NeighbourInterface submodule loaded
before implementations. It must not import NeighbourList, PairPotentials or ACE.
AtomsCalculators properties/units are reused. Do not extend its namespace with
nonexistent functions, use runtime namespace injection, or fork that dependency.

Protocol functions may describe capabilities/preparation requirements. Their
implementations and all execution/allocation belong outside the protocol module.
No mandatory abstract supertype or concrete storage layout for third parties.
Use accessors so implementations are not required to reproduce field names.

## 2. Scope

IN: typed requests; pair and complete-site callbacks; batch callback and CPU loop
fallback; explicit parameters; unit declarations; reusable scratch; current
colored reduction; native adapters; minimal independence/conformance examples;
compatibility and numerical/AD regressions.

OUT: GPU kernels or dependencies, MACE, optimized batching, alternative force
reduction, inverse-incidence/global gradient storage, optimized AD, skin reuse,
general multiple-cutoff composition, new site-energy API, AtomsCalculators proposal.

Keep backend/layout boundaries suitable for a later GPU plan, without building a
general GPU framework now.

## 3. Protocol shape

Names are illustrative; choose a consistent small vocabulary and document final
signatures. Required operations have this meaning:

    evaluate_pair!(out, potential, block, parameters, scratch, properties)
    evaluate_site!(out, potential, site, parameters, scratch, properties)
    evaluate_sites!(out, potential, batch, parameters, scratch, properties)

Mutating callbacks overwrite requested owned slots and return nothing. They do
not reduce into global atom arrays. Requested output storage is supplied by the
backend and can be reused after its contributions are consumed.

Declare minimal queries/hooks for:
- Pair versus complete-site evaluation support.
- Physical search cutoff and required species/layout information.
- Potential-specific prepared data and scratch construction/resizing.

Resolve layout/capability support once during preparation, outside hot loops.
A simple documented method/dispatch contract suffices; avoid a large trait system.
A potential's mathematically required species channels must be distinguishable
from the backend's chosen packing. Never treat arbitrary provider species IDs as
ACE model channel IDs without an explicit map.

Scratch allocation/preparation methods are implemented by the potential package
and called by the backend. They may accept concrete backend/layout/precision and
capacity information. No ACE type checks inside generic traversal.

The CPU implementation may wrap a potential in a concretely typed prepared
evaluator containing model coefficients, angular backend and options. Caller
parameters remain authoritative; refresh derived coefficients when evaluating
updated parameters. Do not expose the whole CPUExecutor to the math callback.

## 4. Properties, parameters and state

Use instances of AtomsCalculators.Energy, Forces and Virial in tuples.
Normalize order and duplicates once at the public boundary. Callbacks specialize
on concrete request tuple types; do not inspect Symbol collections per edge.

Support every nonempty subset of these three properties through the new path.
Forces-only and virial-only require displacement gradients, but need not compute
or return unrequested quantities. If an old callback necessarily calculates
extras, its compatibility adapter may discard them; document that overhead.
An empty public request preserves existing behavior without evaluating sites.
Reject unsupported properties explicitly.

No replacement Energy/Forces/Virial types. Preserve legacy EnergyOnly and
EnergyForcesVirial requests as compatibility entry points mapping to new requests.

Parameters are explicit; existing bundled-parameter APIs adapt at the boundary.
No trainable cutoff or differentiated topology in this phase.
Fixed model structure lives in the potential; caches live in execution/preparation
storage. New callbacks require no model-state argument. Existing AtomsCalculators
and Lux state signatures remain compatible and pass state through as before.

## 5. Geometry, ownership and units

For lattice vectors in columns:
    d_ij = X[:,j] - X[:,i] + H*s_ij
    g_ij = partial E / partial d_ij

Pair callback: energy and g per interaction lane. Backend owns each physical
pair/image once in the native half traversal; no implicit factor of one half.
Legacy fx/fy/fz currently denote the center contribution: verify and adapt their
sign instead of changing them. Existing directed providers can be adapted to half
ownership, or use explicitly documented weights; never silently double energy.

Site callback: E_i and g_ij = partial E_i / partial d_ij per directed edge.
Each selected site contributes once. Its neighborhood is complete within its
cutoff regardless of pair ownership or which centers are selected.
An ACE internal pair term retains its existing site-energy normalization.

Backend reduction:
    F[:,i] += g_ij
    F[:,j] -= g_ij
    W[a,b] -= d_ij[a]*g_ij[b]
Preserve existing virial sign, tensor orientation and public units.
Periodic self images can have zero net atom force but nonzero energy/virial.
Do not deduplicate repeated neighbor indices having different image shifts.

Reuse AC.length_unit(potential) and AC.energy_unit(potential), implementing
declarations on suitable locally owned potential/prepared adapter types and
delegating calculator declarations where appropriate. Numerical input lengths
and cutoffs use the declared length unit; E uses the energy unit; gradients use
energy/length; virial uses energy. No new unit descriptor or hidden default.

Initial native support may remain restricted to Å/eV: explicitly reject other
valid conventions at preparation rather than silently misinterpret them. General
unit conversion across arbitrary parameter structures is not required.
Positions in other convertible Unitful input units must still convert correctly.
Parameters' numeric conventions are the potential's responsibility.

Keep Unitful outside hot loops. Preserve original-order unitful public results
and raw-array boundary conventions. Displacements are computed from current
numeric positions, not stale/detached cached values in a differentiable path.

## 6. Input views and buffer lifetime

Pair views expose species, displacement components, structural/active masks and
lane semantics. Site views expose central identity/species, complete neighbor
species/geometry and masks; species offsets when required by the prepared layout.
The backend retains destination identities/image labels needed for reduction.
Distance values may be supplied additionally; vectors are mandatory for ACE.

Batch views expose site count and access to each complete environment and output;
variable neighborhood sizes use offsets/ranges, with padding explicitly masked.
Do not mandate Matrix, Vector-of-Vectors, SIMD.Vec, a fixed device type or fields
in the protocol. Native views can expose SIMD efficiently through specialized
accessors. IDs and any permutations retain an explicit map to original atoms.

Views are borrowed until the callback finishes, or until an immediately ensuing
backend consumption step, as documented. A potential must not retain them.
Requested active outputs are overwritten; inactive gradient outputs are zeroed
or provably never consumed. Apply safe dummy geometry BEFORE evaluating functions
that could otherwise produce NaN/division by zero. Masking after NaN is inadequate.

## 7. Execution and batching

Backend types own preparation and execute!/evaluate! methods; their algorithms
do not belong in the protocol. Retain the current public evaluator/calculator
behavior through adapters. Public evaluations may rebuild topology every time.

batch_size is a positive configurable maximum independent of SIMD width, thread
count and scratch-slot count. batch_size=1 calls the direct single-site path.
For larger batches, a CPU implementation/utility fallback loops over site
callbacks; the pure protocol does not implement this fallback.

Batches contain only selected centers and retain full neighborhoods. Serial
execution follows the documented selected-center order. Parallel direct writes
respect existing complete endpoint footprints and phase barriers. Batch within
compatible colored work, retain exclusive slot-owned scratch and handle a smaller
tail. If batch output buffers are used, size them at preparation and consume all
site outputs before reuse. Never alias multiple unconsumed outputs.

Retain one shared final force array, no full force array per worker.
Use only bounded per-slot/site/batch scratch for the current colored path.
Do not build a global edge-gradient array or add gather reduction now.
Keep energy/virial partial ownership race-free.

## 8. Implementation mapping

Current ACE/environments.jl mixes channel mapping/cutoff selection with full
connectivity, grouping and coloring. Separate mathematical preparation from
generic storage/scheduling while preserving FullEnvironment compatibility.

Current ACE/cpu.jl::cpu_site! mixes mathematics and global accumulation.
Extract the math into a site callback; write gradients into local supplied
buffers, then reduce in NeighbourList. Retain radial/angular/projection kernels
and reference paths. ACE-specific SiteScratch stays in ACE.

PairPotentials/evaluation.jl owns current traversal/reduction. Move/delegate its
generic execution responsibilities to NeighbourList as needed. Keep pair model
binding and adapters in PairPotentials. Preserve fused block traversal rather
than forcing materialized atom lists.

Keep the root thin and dependencies one-way:
NeighbourInterface <- NeighbourList
NeighbourInterface <- PairPotentials / ACE
Potential models and calculators may use NeighbourList, not the reverse.
Avoid circular includes and concrete ACE dispatch checks in neighbor loops.

## 9. AD and compatibility

Preserve supported energy->positions/parameters, forces->positions/parameters
and virial->parameters derivatives. Preserve documented limitations such as
unsupported virial->positions. Retain current serial/ForwardDiff/ChainRules
reference path; optimizing AD or adding derivative rules for every new mutating
callback is outside scope. Do not claim the new prepared executor is automatically
differentiable.

Calculator/Lux differentiation must still reach the supported functional entry
points. Preserve input snapshots, delayed/repeated pullback correctness,
workspace isolation and parameter authority. No new precision casts that erase
dual numbers. Preserve custom radial functions/source invocation safeguards,
projection/masking/pruning, selected centers, fitting and persistence.

Persist models/parameters as before, not scratch/executors/views. If environment
types move, maintain old qualified names through aliases/adapters where practical.
No package rename, UUID/author/license changes or unrelated cleanup.

