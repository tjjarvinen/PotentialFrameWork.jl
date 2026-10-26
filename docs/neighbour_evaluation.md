# Experimental CPU neighborhood evaluation

`PotentialFrameWork.NeighbourInterface` (NI) is a local experiment using
AtomsCalculators (AC) property instances and unit declarations. AtomsCalculators
is unchanged. The protocol imposes neither a potential supertype nor view fields.
A runnable third-party pair and angular-site implementation is in
[`examples/neighbour_evaluation.jl`](../examples/neighbour_evaluation.jl).

Implement these methods on your own potential/prepared types:

```julia
NI.evaluation_kind(p)                         # :pair or :site
NI.search_cutoff(p)                           # numeric physical radius
NI.species_channels(p)                        # ordered species, or nothing
AC.length_unit(p); AC.energy_unit(p)          # explicit declarations
NI.prepare_potential(p, species; precision, width, kwargs...)  # concrete evaluator
NI.make_scratch(prepared, capacity, precision)                # exclusive scratch
NI.evaluate_pair!(out, prepared, block, parameters, scratch, properties)
NI.evaluate_site!(out, prepared, site, parameters, scratch, properties)
NI.evaluate_sites!(out, prepared, batch, parameters, scratch, properties) # optional
NI.update_prepared!(prepared, parameters)     # optional derived coefficient refresh
```

Implement the callback matching `evaluation_kind`. Parameters are authoritative
on every execution. Callbacks have no required model-state argument. Model state
continues through existing calculator/Lux adapters. Structural cutoffs are fixed;
trainable topology and skin reuse are outside this interface. Preparation resolves
capabilities/species/layout once; unsupported dispatch/layouts must reject there.
Native execution accepts declared Å/eV only and rejects other conventions.
Convertible Unitful system positions/cells are converted at extraction. Numeric
parameters follow the potential's declared convention; no arbitrary tree conversion
is attempted.

Access views through these functions:

| View | Accessors |
| --- | --- |
| Pair | **`distance(block)`** (scalar/SIMD distance; primary radial input), `central_species(block)`, `neighbour_species(block)` (lane tuple), `neighbour_count(block)`, optional `displacement(block)` (x/y/z scalar or SIMD components), `active(block)`, `structural(block)`; indexed species/masks are also available |
| Site | `central_index(site)` (original atom), `central_species(site)`, `neighbour_count(site)`, `neighbour_species(site,k)`, `displacement(site,k)` (three-vector), `active(site,k)`, `structural(site,k)`, `species_range(site,z)` (local padded range for declared channel z) |
| Batch | `site_count(batch)`, `site_at(batch,k)`, `output_at(out,k)` |
| Output | `set_energy!(out,E)`, **`set_radial_derivative!(out,dEdr)`** (usual pair output), optional `set_gradient!(out,g)` (Cartesian pair component tuple), `set_gradient!(out,k,g)` (site edge), `output_energy(out)`, `output_radial_derivative(out)` (radial pair), `output_gradient(out[,k])` (Cartesian pair/site) |

Site channel IDs follow the potential's declared ordering, through an explicit
identity map from provider species. They are not arbitrary provider IDs. ACE
retains its strict element-only conversion and rejects named/isotopic channels.
With `species_channels(p) === nothing`, native packing uses distinct provider
species in encounter order; callbacks requiring mathematical channel ordering
must declare it. Padding is inactive and has safe geometry **before** evaluation.
Do not consume inactive outputs; their contents are unspecified.

The main pair implementation is radial: obtain `r = NI.distance(block)` and
write `NI.set_radial_derivative!(out, dEdr)` for force/virial requests. Callbacks
need no orientation. The backend owns the displacement and converts `dEdr` to
`g = (dEdr/r)*d` before force and virial reduction. Distances are computed from
current safe block geometry, without requiring additional stored arrays or square
roots for legacy Cartesian callbacks. Native `RadialPotential` and `SmoothPair`
use this path. `NI.displacement(block)` and Cartesian `NI.set_gradient!` remain
optional for orientation-dependent pairs. Set one derivative representation per
callback; setters select which representation the backend consumes. At physical
r=0, radial models must have a regular origin (`dEdr=0`) or explicitly reject the
geometry. Backend division uses a safe denominator, preserving the established
regular-origin convention. The historical callable pair adapters remain available
for the functional AD and legacy entry points.

Callbacks overwrite requested outputs and return `nothing`. They may not reduce
into global atom arrays or retain borrowed views. Output slots remain valid through
immediate backend consumption; every pending batch site has distinct storage.
Energy is per owned physical pair/image or per selected complete center, with no
implicit half factor. The optional Cartesian gradient is `g = ∂E/∂(Xj-Xi+H*s)`;
the default radial output is `dE/dr`, with the same sign convention. Backend reduction is
`Fi += g`, `Fj -= g`, `W[a,b] -= d[a]*g[b]`. Different periodic images are never
merged by neighbor index. Nonzero self-images can have zero net force and nonzero
energy/virial. Legacy `fx/fy/fz` already have this gradient sign. Legacy callbacks lacking unit
methods retain their established Å/eV contract through the adapter; new protocol
potentials must declare units explicitly.

`NeighbourList` owns full topology/packing, views, buffers, loop fallback, colored
scheduling and reduction. `ACE.FullEnvironment` remains an alias for the moved
type. ACE owns channel/cutoff policy, scratch, kernels and derived coefficients.
PairPotentials owns legacy callback/model binding. No ACE checks occur in generic
backend traversal. The protocol module contains declarations and contracts only.

## Native entry points and batching

For pair callbacks, use the existing `evaluate(sys,p; request=(AC.Energy(),...),
parameters=ps, width=4)` / `evaluate!(workspace,sys,p; ...)` boundary. Legacy
`EnergyOnly()` and `EnergyForcesVirial()` remain supported. Existing bundled pair
callbacks and scalar differentiable adapters retain their APIs.

For a protocol site potential:

```julia
input = NeighbourList.prepare_site_execution(p, sys;
    precision=Float64, width=4, slots=Threads.nthreads(), batch_size=2)
NeighbourList.execute!(input.execution, input.positions, input.cell, ps;
    properties=(AC.Energy(), AC.Forces(), AC.Virial()), centers=nothing,
    parallel=false)
y = NeighbourList.site_result(input.execution) # owning original-order total E/F/W
```

This is an executor interface returning totals, not a public site-energy API.
Use current original-order numeric positions, and column-vector cell H with open
generators zeroed. Prepared topology must remain valid; rebuild preparation after
connectivity/species changes. Environment version/capacity changes reject stale
execution. No automatic skin or topology update is provided by the numeric call.

`ACE.CPUExecutor` and `ACE.ACECalculator` also accept `batch_size` at construction;
`ACE.execute_prepared!` accepts `properties`. Default batch size is 1. It calls the
site callback directly. Larger positive integer sizes call `evaluate_sites!` with
a smaller tail; the CPU fallback loops over complete sites using one scratch slot
and distinct output buffers. A potential-specific batch method specializes on its
prepared type without requiring native batch fields. The default no-op
`update_prepared!` and batch loop are implementation methods in NeighbourList,
not protocol algorithms. They also work with other accessor-based views.

Serial evaluation follows the requested original center order. Parallel execution
batches within colored phases and assigns scratch by exclusive task slot, never
thread ID. All complete center/neighbor force destinations are colored. There is
one shared sorted force array; bounded per-slot gradients cost
`3*sizeof(T)*capacity*batch_size*slots` bytes, plus small output/center metadata.
No global edge-gradient array or full force array per worker is introduced.
Callbacks add a gradient-buffer write/read and backend displacement reconstruction
relative to the formerly fused ACE reduction. Batch fallback is a correctness
feature; no batching throughput gain is claimed.

Every nonempty subset of AC.Energy/Forces/Virial is supported. Boundary requests
normalize order/duplicates once for typed callbacks; unsupported requests reject.
Empty requests evaluate no interactions. Existing calculator property-result key
order, state, units, additive destinations and owning outputs remain compatible.
Legacy pair adapters may compute/discard extra quantities for gradient requests;
ACE currently retains its joint radial/angular derivative kernels for energy-only
prepared requests. Unrequested properties are not accumulated or consumed.

## AD, independence and limitations

The functional `pair_evaluate`/`ace_evaluate` and existing ForwardDiff/ChainRules/
Zygote/Lux adapters retain supported energy→positions/parameters,
forces→positions/parameters and virial→parameters products, including saved input
ownership. Documented unsupported derivatives still reject. New mutating prepared
executors have no new AD rules and must not be assumed differentiable. Existing
persistence saves model/parameters/policies, never executors or views.

`test/neighbour_evaluation.jl` supplies different-storage views and an independent
O(N²) open-cell backend running the same pair and true angular callbacks, with
finite-difference force/strain checks. It explicitly rejects periodic geometry.
Native paths keep complete skew/slab/wire/multi-image support and optional
NeighbourLists query interoperability. No second optimized provider executor,
GPU, MACE, alternative reduction, optimized contractions/AD, or upstream protocol
proposal is part of this experiment.

Additional standalone validation compares all 33 local saved full-precision
ACEpotentials predictions with prepared CPU evaluation across batch sizes 1/2/5
and widths 1/4: `test/neighbour_saved_reference.jl`. It validates source/geometry
hashes and reads the supplied files without loading/refitting the reference model.
Prepared timings/allocation and inference inspection can be reproduced with
`benchmark/neighbour_evaluation.jl` and `benchmark/inspect_neighbour_evaluation.jl`.
