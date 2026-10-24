# CPU ACE inspection and transformation interface

All APIs below live in `PotentialFrameWork.ACE` (`A` in the examples). They extend
resolved imports and native checkpoints. Existing calculator, numerical E/F/W,
Unitful labels, AD ownership and directed periodic-image conventions still apply.
No new site-energy property or general ACE basis generator is introduced.

## Inspect and select complete invariant features

```julia
using PotentialFrameWork, SpheriCart
A = PotentialFrameWork.ACE
m, ps, st = A.load_resolved("data/ace_ref/resolved_model_v1.toml")
summary = A.basis_summary(m, ps)       # owning arrays, no plotting dependency
rows = A.descriptors(m)              # actual coupling/product dependencies
selection = A.BasisSelection(m; indices=[1, 4, 9, 17])
# Alternatives: ids=[rows[1].id, rows[4].id], orders=1:2,
#               max_index_weight=8 (criteria intersect)
masked = A.mask_basis(m, selection)
reduced = A.prune_basis(m, ps, selection)
small_model, small_ps = reduced.model, reduced.parameters
```

A descriptor ID is `(source_basis_namespace, original_row)`. Metadata includes
current/original indices, each coupled product with current/original density and
radial dependencies and exact harmonic (l,m) tuples, source coefficient IDs,
actual correlation order set, coefficient indices by central channel, mask and
species-order policy. Mixed-order rows retain their complete order set. An order
criterion retains a row only if **all** its orders satisfy the criterion. Unknown
custom index weights are `nothing`; selection by those undefined weights fails.
The index weight is algebraic sum(n+l), not a physical polynomial degree.
Nominal body order adds the center to correlation order; repeated atom factors
mean this is not a distinct-atom count. Metadata arrays are owned copies.

Masks preserve dimensions, ordering and unmasked ps. They are applied on every
scalar/CPU/fitting/AD call, including after optimizer or readout installation.
Disabled entries have zero parameter sensitivity. Reapply an all-row selection
to remove a mask. Pruning compacts invariant rows, structural product columns,
density indices, imported radial channels and nonlinear radial scales. It never
uses fitted coefficient zeros to choose structure. Original IDs survive repeated
reductions and persistence. Pair coefficients/scales and references remain unless
explicitly disabled with `A.components(m; pair=false, reference=false)`.

`reduced.maps` maps current invariant/product/density/radial/harmonic dependencies
to the immediate parent. `reduced.coefficient_map` is T in `c_parent=T*c_reduced`;
`A.mapped_prior(R,T)` returns `R*T`, retaining the original prior action. For
`fit_pair=false`, use `A.coefficient_map(m, selection.indices; fit_pair=false)`.
`A.restrict_parameters(ps, reduced.maps)` provides a differentiable numerical
gather for lifting gradients back into parent readout/radial-scale coordinates.
Its unchanged subtrees are shared; `prune_basis` returns owning parameter copies.
Arbitrary coordinate maps are allowed for prior transfer; they need not remain
diagonal. The ACEfit extension explicitly augments rectangular prior operators
with matching zero targets, preserving this objective despite ACEfit 0.3.1
QR assuming square priors. `install_readout(ps,c)` retains its existing API: mathematics and masks
belong to the model, rather than overwriting optimizer parameters.

`A.restrict_basis(m,ps,selection)` has the same reduction semantics.
`A.restrict_basis(path; indices=..., orders=..., max_index_weight=...)` reads an
existing resolved specification and returns the restricted model, parameters,
state and maps **before constructing any executor or numerical workspace**.
It reads the existing full basis specification; it does not generate new ACE
invariants. Angular execution still computes full harmonics through the existing
lmax. Pruned custom radial adapters can retain computation of their parent's
channels; their maps are compact but that callable work is not claimed removed.

## Original-order centers and neighbor-species factor orders

```julia
projected = A.project_species(masked, [:Ti, :Al], 1:2)
input = A.ace_inputs(sys, projected)
r = A.ace_evaluate(input.positions, projected, ps,
                   input.environment, input.cell; centers=[3, 1])
design = A.design_matrix(input.positions, projected, ps,
                          input.environment, input.cell; centers=[3, 1])
calc = A.ACECalculator(projected, ps, st)
# AtomsCalculators.energy_forces_virial(sys, calc; centers=[3, 1])
```

`centers=nothing` means all original atom indices. Ranges, ordered unique vectors
and empty selections are supported; duplicates and invalid indices fail.
Selections never slice the system: every center uses its full environment, and
forces accumulate onto **all** affected atoms. The output force matrix/vector is
always full-system original order. Empty selections return zero E/F/W. Imported
pair terms use their existing directed pair-site partition, not complete
unordered pairs. CPU phases filter center tasks while retaining full write
footprints and one shared sorted force array. Preparation still builds the full
system environment; this has not become a region-only neighbor builder.

One species/group has one combined neighbor-factor counter. The center's species
never contributes. For each density projection A=U+V (selected species in V),
products are expanded as coefficients of `prod(U[a]+t*V[a])`. Orders in the given
inclusive range are summed with the **same original readout coefficient**.
Production uses truncated convolution and a division-free analytic reverse, with
preallocated CPU polynomial/adjoint scratch. Intermediate lower orders are kept
even when the minimum requested order is high. There is no exponential chemical
basis or independent coefficient per k. Complete coupled invariants, rotation
symmetry and permutation symmetry are preserved.

An empty species group has V=0: order zero reproduces the original model and
positive orders vanish. Unknown/named/isotopic species follow the existing strict
element-channel rejection policy. Malformed/negative ranges fail; orders above
the product order vanish. `0:max_order` recovers the original model. References
have order zero. Pair edges have order one for selected neighbor species and zero
otherwise; component flags also apply. Thus positive-order models remove the
reference baseline. `site_diagnostics` remains available for individual centers;
no new per-site-energy API is added. Projections are algebraic decompositions,
not isolated-subsystem or QM/MM physics.

Centers are explicit per-call discrete metadata. For Lux, add `centers` to the
input NamedTuple; models contain masks/species policies. Numerical AD keeps
centers/topology fixed, owns saved inputs and accumulates all repeated destinations.
E->X/ps, F->X/ps, W->ps and the existing W->X products remain supported with
ForwardDiff and Zygote. Arbitrary reverse derivative order and cell AD are not
claimed. Fitting uses identical selected/transformed rows, baselines and masks;
Unitful labels are still converted to eV/eV/Å before baseline subtraction.

## Custom radial callables

`RadialMetadata(spec; cutoff, index_weights=nothing, simd=false, ad=false,
origin=:error, envelope=:included, channel_dependence=:central_neighbor)` declares
output size/order via (n,l) labels, fixed Å cutoff, algebraic weights, capabilities
and origin/envelope/channel policies. The callable owns its cutoff envelope
exactly once. Radial/angular correspondence remains that of the existing basis;
`replace_radial` rejects ordering or species-matrix shape mismatches.

The in-place protocol is
`f(values, derivatives_or_nothing, distance, custom_ps, config)`; functions and
callable structs are accepted. `A.CallableRadial(f, metadata; config=...)` stores
structure/configuration separately from ps. `A.radial_evaluate!` writes caller
buffers. Values-only calls pass `nothing` as derivative storage. `radial_ed` and
`radial_values` are allocating scalar conveniences. Imported spline channels keep
their specialized SIMD Horner implementation and support the same scalar protocol.

```julia
meta = A.RadialMetadata(m.radial_spec; cutoff=5.5,
    index_weights=sum.(m.radial_spec), simd=true, ad=true, origin=:regular)
config = (powers=[n-1 for (n,l) in m.radial_spec], cutoff=5.5)
ch = A.CallableRadial(A.SmoothExponential(), meta; config,
    persistence=(route="registered", id="PotentialFrameWork.smooth_exponential.v1"))
custom = A.replace_radial(m, ps, fill(ch, 2, 2);
    parameters=(rate=[0.7], amplitude=[0.2]))
# custom.parameters.custom_radial is the authoritative custom parameter subtree.
# pair=true uses the pair specification and parameters.custom_pair instead.
```

The portable exponential/power example has trainable rate/amplitude and a smooth
quartic cutoff. `examples/ace_radial_source.jl` contains a self-contained generic
implementation; `examples/ace_interface.jl` demonstrates source-based use.
Species-specific configurations use an explicit homogeneous matrix of channels;
central/neighbor channel selection occurs outside the distance/channel kernel.
The callable must produce the declared order for both scalar and SIMD.Vec inputs.
Each Vec lane is an independent distance; dynamic integer powers in SIMD.jl may
need explicit multiplication. Unsupported arbitrary SIMD AD is rejected; CPU
callable execution currently uses scalar lanes, while direct generic Vec calls
can be validated separately. `SubsetRadial` pruning is a correctness adapter
that evaluates parent callable channels before gathering; it may allocate.

`CallableRadial(values_only!, meta; derivative=:ad)` uses the optional ForwardDiff
scalar distance-derivative adapter. This fallback allocates channel temporaries
and reevaluates values per output derivative. `A.allocating_radial(f,meta; ...)`
adapts an allocating `(values, derivatives)` function; `derivative=:ad` instead
adapts an allocating values-only function. Both are explicitly convenience paths.
Callable pullback snapshots copy mutable configuration and retain immutable code
references (including captured module bindings). Noncopyable mutable captures
must move into explicit config/ps rather than borrowing caller state.
Joint analytic derivatives must remain generic for second distance/mixed parameter
products needed by force losses. A derivative declaration alone is not evidence
of generic AD support: independent force-loss checks are required. Undeclared
custom AD is rejected by the reverse-product boundary.

`A.radial_capabilities(ch, custom_ps)` reports scalar/Vec value-and-first-derivative
checks, an independent `training_ad_declared` field and diagnostic reasons. With
ForwardDiff loaded, `training_ad_checked`/`training_ad_local` also probe local
first-derivative consistency and second/mixed derivative execution at samples;
these are separate from end-to-end model force-loss validation. The helper
never reports verified GPU support. Plot samples can be obtained without plotting
from `A.radial_samples(m; channels=1:4, parameters=custom_ps)`. Unknown custom index
weights require an explicit prior or `smoothness_prior(m; unknown=:neutral)`, an
explicit neutral identity prior, rather than an invented smoothness degree.

The optional `examples/acepotentials_radial_adapter.jl` checks actual pinned
ACEpotentials 0.10.2 `Models.evaluate_ed` and its normalization/order/envelope.
Run it with `julia --project=data/ace_ref examples/acepotentials_radial_adapter.jl`.
It loads the saved model without fitting. This scalar allocating adapter does not
claim SIMD/AD/persistence capabilities or compatibility with arbitrary full models.

## Checkpoints and explicit approximation

`save_calculator` writes version 2, while old native v1 and resolved v1 files still
load. It saves authoritative ps/st, complete numerical basis, masks, provenance,
component and species-order policies. Environments, scratch and precontractions
are rebuilt. Per-call center choices are not mutable calculator state.

Custom radial persistence routes:

1. Register a trusted versioned implementation with
   `register_radial!(id, (config,metadata)->evaluator)`. The checkpoint stores ID,
   config and numerical references; load its package/extension first. The built-in
   exponential above has a shipped registration.
2. `source_radial(source_text, "factory", metadata; config, allow_code=true)`
   calls a factory `(config,metadata)->evaluator`. Read source files explicitly
   with `read(path,String)`; include all definitions rather than captured globals.
   Bundles store source/hash/interface version, configuration, optional exact
   dependency versions and value/derivative references. Loading any checkpoint
   containing source defaults to refusal before source parsing/evaluation:
   `ACECalculator(path; allow_code=true)` is explicit trust. Source runs in a
   fresh namespace; this is **not** a sandbox. Hash/reference checks establish
   identity/numerical consistency, not safety. No packages are installed and no
   saved shell commands are executed by the loader. Public execution boundaries
   handle world age; radial loops contain no `invokelatest` calls.
3. `tabulate_radial(ch, custom_ps; domain=(rmin, cutoff), nodes=1001,
   atol=1e-6, derivative_atol=1e-5)` is an **explicit approximation**. It returns
   `.channel` and `.report`, using raw-distance cubic Hermite interpolation and
   independent off-grid checks. Replace channels explicitly with this table and
   empty custom parameters; persist it without custom code. The lower domain
   bound is enforced; beyond the upper cutoff outputs are zero. Tables are C1,
   generally not C2 at knots: vibrations/force-loss second derivatives can jump.
   Original custom-parameter gradients do not survive this fixed snapshot.

Unsupported callable saves fail before replacing the destination. Source loading
checks dependency versions, hash, factory, output size, scalar/declared SIMD
capabilities and numerical references (atol 5e-10, rtol 5e-9). These are checks at
finite samples, not proofs for all inputs. Custom source constructors execute
only when explicitly permitted. Tabulation reports measured value/derivative
errors and requested tolerances; no equivalence beyond them is claimed.

## Optional terminal plots and backend boundaries

```julia
using UnicodePlots
A.plot_basis(m, ps; kind=:orders)
A.plot_basis(m, ps; kind=:weights)
A.plot_basis(m, ps; kind=:coefficients, indices=1:30)
A.plot_basis(m, ps; kind=:coupling, indices=1:30)
A.plot_basis(m, ps; kind=:counts, indices=1:30)
A.plot_basis(m; kind=:radial, channels=1:4)
```

Each call returns a UnicodePlots object; it never displays implicitly. Empty bases
are labeled with zero counts. Coupling/coefficients and radial curves default to
bounded selections. Coefficient magnitudes depend on normalization and coordinate
representation; they are exploratory, not an automatic importance/pruning rule.

| Path | Requirements / status |
| --- | --- |
| CPU scalar | Generic numeric values and joint distance derivatives; explicit ps/config; origin/cutoff policy. |
| CPU imported SIMD | Existing spline/Horner fast path and safe inactive-lane handling retained. |
| CPU callable lanes | Exclusive scratch and in-place analytic callable buffers; subsets/AD convenience adapters may allocate. |
| Direct callable Vec | Generic independent-lane arithmetic, declared/tested separately; no arbitrary Vec AD promise. |
| Training AD | Generic first derivatives with second/mixed products; tested ForwardDiff/Zygote boundaries; topology/centers discrete. |
| Future GPU | Deferred/unverified. Requires device-compilable scalar operations, concrete device-compatible state, no host-only calls/unsupported allocations and device compilation tests. Vec acceptance is not GPU support. |

Inactive CPU lanes are skipped or given safe dummy directions before transforms;
zero times NaN is not protection. Coincident physical atoms retain the existing
angular/origin error policy. Public results own storage; prepared scratch is
borrowed only inside exclusive execution slots. No force array per worker is
introduced. Model sizes, callable-adapter costs and projection overhead must be
measured separately from full environment preparation and owning output costs.
