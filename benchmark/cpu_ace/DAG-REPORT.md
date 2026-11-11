# CPU ACE product DAG

10 November 2026. This report measures the incremental step from division-free
prefix/reverse products to a deterministic shared product DAG. The stage-one
report and all its original results remain unchanged. Raw measurements and
[complete tables](results/2026-11-10-dag/OBSERVATIONS.md) accompany the
[reproduction instructions](results/2026-11-10-dag/REPRODUCE.md).

## Measured result

On AMD Ryzen 9 9950X, Julia 1.13.3, one pinned physical core and BLAS=1,
large TiAl at 1,024 atoms improves **139.909 → 68.453 ms prepared (2.04×)**
and **155.992 → 83.257 ms public (1.87×)** versus prefix/reverse. Polynomial
reference is 594.411 ms prepared: the incremental DAG gain is distinct from the
previous polynomial-to-prefix improvement. At 11,664 atoms, Float64 prepared
improves **1,663.012 → 775.673 ms (2.14×)** and public **1,845.275 → 975.139 ms
(1.89×)**. Float32 gains there are 2.16× prepared / 2.05× public.

At 1,024 atoms, published Ti and Ti–Al–V Float64 prepared gains are 1.58× / 2.16×.
Small TiAl and medium Ti gain about 5–7% prepared and 3–6% public in the serial
grid. Small TiAl's 11,664-atom Float64 public median is effectively unchanged:
313.353 → 313.349 ms. Float32 is not universally faster; the published Ti–Al–V
polynomial Float32 prepared reference is slower than its Float64 reference.

All **612 BenchmarkTools trials reached ten warmed samples**, with no
under-sampled cases. Every observation, actual count, min/max spread, GC time and
host allocation volume is retained in TOML; `summary.json` audits the trials.
This is one controlled benchmark sequence, with independent one-core repeats in
the scaling job and an immutable-source repeat, not a statistical population of
hosts. Expensive trials ran sequentially. Compilation/correctness warmup preceded
production samples; the helper's `warmup_s` may already be compiled and is not
presented as compilation latency.

The first-site large TiAl Float64 product/reverse kernel improves 104.471 →
32.845 µs (3.18×), comprising 9.785 µs forward and 23.110 µs reverse when measured
separately. Its unchanged readout is 0.870 µs. Small-model products improve 2.41×,
but their complete prepared gain is only 1.06× at 128 atoms. The complete tables
include every model/precision/size and both kernel and workload boundaries.

Large TiAl's simple plan has 165 leaves, 13,919 original/unique terminals,
16,229 multiply nodes, 2,334 intermediate-only nodes and depth 3. Its independent
nontrivial-product sharing ratio is 2.23. Forward multiplication estimates fall
from 50,080 to 16,229 and reverse from 100,160 to 32,458. The small and medium
models have sharing ratios 1.48 / 1.40 and much smaller product fractions. Product
fraction, sharing, seed/clear costs and the rest of the force stack together
explain why sharing ratio alone does not predict complete-call gains.

## Construction, memory and regressions

The large TiAl Float64 plan builds in 1.149 ms, allocating 7,333,584 bytes during
construction; retained plan-array payload is 771,656 bytes. Its per-slot DAG
numerical payload is 262,304 bytes, including leaves; retained polynomial tables
are 400 bytes, and total numerical slot scratch at 128/1,024 atoms is 484,368
bytes. Float32 halves numerical payloads, including the 131,152-byte DAG scratch,
while retaining identical integer plan arrays. Published Ti–Al–V retains
2,534,640 plan-array bytes and 921,616 DAG numerical bytes per Float64 slot;
its measured build is 4.459 ms / 25,860,176 allocated bytes. Full counts for all
fixtures, including fallback and total scratch, are in the observation tables.

At 128 atoms the isolated plan build is less than one saved prepared-call time
for every measured model/precision. This is a plan-only amortization estimate;
public-call savings can be negligible, as in large-system small TiAl Float64.
Ordinary calculator reuse pays graph construction once. Fresh-calculator calls
also include construction/scratch/owning outputs: large TiAl Float64 at 1,024
atoms improves 157.932 → 84.434 ms on the integrated checkout, and at 11,664
1,850.498 → 975.889 ms.

The immutable stage-one source exposes costs common to the new switches. Its
large TiAl Float64 prepared median at 1,024 atoms is 140.958 ms versus 139.909 ms
for the integrated prefix path; public medians are 156.013 / 155.992 ms. Fresh
prefix calculator calls **regress 155.925 → 157.932 ms** at 1,024 atoms, and
20.091 → 20.914 ms at 128 atoms. These measured construction regressions are not
hidden by same-checkout comparisons. The new plan/scratch exists for all switches
to support changing algorithms without graph construction during execution.

Every serial numerical prepared trial in the main and 11,664-atom grids allocates
**zero warmed host bytes**. `parallel=true` allocates task/scheduling storage even
with one core: large TiAl allocates 311,344 / 460,768 / 1,005,936 bytes at 1/4/16
cores. At sixteen slots its numerical scratch payload is 7,749,888 bytes; the one
shared force payload remains 24,576 bytes for 1,024 atoms. These are separate from
plan payload, construction allocation volume and allocator peaks.

## Scaling and remaining bottlenecks

Large TiAl 1,024-atom Float64 prefix → DAG prepared medians at 1/4/16 physical
cores are 139.537 → 68.523, 55.215 → 35.911, and 38.484 → 21.065 ms. Incremental
gains are 2.04× / 1.54× / 1.83×; public gains are 1.86× / 1.35× / 1.54×. DAG
prepared scaling from one to sixteen cores is 3.25×, not linear scaling.

Small TiAl public medians slightly regress at four cores (21.397 → 21.485 ms)
and sixteen (19.087 → 19.100 ms). These differences are small relative to observed
variation, so they do not establish a statistically significant loss. Parallel
small-model prepared timings are particularly variable: sixteen-core DAG samples
span 2.063–6.610 ms and prefix 4.187–8.729 ms. Do not interpret the median's 1.31×
ratio as a stable universal gain. Affinity is controlled and BLAS is one thread;
no scheduler or coloring changes were made to improve these results.

Undisturbed serial workload profiles anchor the next-bottleneck assessment.
In the 1,024-atom large TiAl Float64 DAG prepared profile, the product call has
1,058 of 2,008 snapshots (~53%), with its reverse call at 771 (~38% of total).
Array reads/writes and bounds checks remain visible beneath reverse propagation.
For small TiAl Float64 public evaluation, preparation accounts for 449 of 756
snapshots (~59%). These are sampled call-stack observations from one profiling
run per workload, not additive instrumented-stage percentages or cache counters.

Concrete production forward/reverse and wrapper LLVM/native/warntype dumps were
inspected in Float32 and Float64, beyond the wrapper's two call instructions.
Float32 callees use float arithmetic, with no float-to-double extension or double
arithmetic observed; Float64 uses double arithmetic. No direct GC allocator calls
were found in these bodies, and all six warmed callee allocation probes report
zero bytes. Index checks/error branches remain (11/12 static bounds-error call
sites in forward/reverse); error helpers can allocate if an invalid index throws.
Static counts are not runtime cost percentages. `perf` is not installed, as the
saved probe records; no hardware cache/bandwidth claim follows from these data.

## Implemented boundary

`:dag` is explicitly selectable through CPUExecutor, ACECalculator, execution
keywords and generic NeighbourInterface preparation. Polynomial `:reference`,
independent prefix/reverse `:optimized`, and restricted-order polynomial fallback
remain runnable. Auto remains prefix/reverse. Functional energy/force/virial
parameter/position differentiation and equivariant feature execution retain
polynomial projection and their existing derivative rules. The new mutation
path supplies inference density adjoints; it adds no training AD capability.

The plan sorts each scalar product's factors, preserving multiplicity, and interns
operand pairs while walking prefixes. Leaves and multiplication nodes occupy
disjoint index ranges; internal nodes are topologically ordered. Every original
coordinate maps to a terminal, including duplicate products and terminals used
by other products. All weights seed terminals additively, then reverse traversal
adds both operand contributions, including equal operands. There is no division,
zero-weight pruning, parameter specialization, coupling/coordinate modification,
fast-math, contraction rewrite, or bounds-check optimization.

One shared read-only structural plan belongs to preparation. Each exclusive task
slot owns node values and adjoints. Reverse clears all node adjoints every call;
forward overwrites all node values and outputs. Authoritative coefficients/masks
refresh existing product weights on every execution. Structural snapshots detect
changed product indices/density count at the executor/update boundary; changed
models need new preparation. Plan construction never occurs in a production site
loop. Existing PreparedACE and array scratch constructors remain compatible.
Legacy array scratch may resize its numerical buffers on its first DAG call;
normal executor and generic preparation allocate them before site execution.

The executor retains polynomial fallback scratch and its one sorted shared force
array. Graph arrays, numerical scratch, observed host allocations, construction
allocations and allocator peaks are different quantities. This implementation
allocates the graph and DAG scratch during preparation even when prefix/reference
is selected, so immutable-source comparisons are necessary to expose shared
integration overhead. No transient graph/scratch is added to model persistence.

## Measurement boundaries

Fixed width 4, SIMD radial enabled, fusion disabled, SpheriCart, BLAS=1; both
precisions use identical product settings. Perturbed nonsymmetric bcc geometry,
models and cached independent oracles are unchanged. Product diagnostics use
actual first-site U/V state and the real production callees. DAG forward and
reverse are measured separately and together; prefix's fused per-product
forward/reverse remains one combined diagnostic. Readout is measured separately;
production energy accumulation remains unchanged. These isolated measurements
perturb caches and are not additive whole-call runtime percentages.

Prepared E/F/V includes weight refresh, topology checking, products, gradients,
scheduling and shared reduction, with prepared geometry metadata. Public calls
reuse the calculator and graph, but perform boundary extraction, fresh neighbor
construction/packing/coloring, executor preparation and owning unitful outputs.
`construct_and_public` also constructs a fresh calculator/graph/scratch each
sample. Plan build and executor construction are measured separately. There is
no hidden mutable global cache or CPU skin reuse.

Operation estimates count prefix's multiplication by the initial identity:
forward one and reverse two multiplies per original factor. DAG forward uses one
and reverse two per multiplication node. Prefix forms U+V twice per factor and
adds two density adjoints per factor; DAG forms U+V once per leaf, clears all
node adjoints, accumulates one seed per original product, adds two adjoints per
multiplication node, then copies leaf adjoints to both density outputs. The
reported sharing ratio divides independent nontrivial-product multiplies
(sum(order-1)) by DAG node count; it is not an expected speedup.
Plan payload includes left/right/terminal arrays and copied original structural
indices used for compatibility checks. It excludes object headers and allocator
metadata. Per-slot DAG scratch includes separate leaf values and adjoints, and
total scratch includes radial/angular/neighborhood arrays and retained polynomial
fallback. All-slot scratch and shared forces are recorded for each thread count.
Neither payload nor measured allocation volume is an allocator peak.

## Validation

Full `test/runtests.jl` passes **56,658 assertions in 62 testsets per process** on
Julia **1.13.3 and 1.10.12**, four threads; 2,065 assertions are added relative to
the recorded stage-one suite. Focused product/integration tests pass 2,700
assertions. Separate ACEfit 0.3.1 validation passes four assertions. CPU-only
loading passes with AMDGPU/CUDA/SpheriCart/UnicodePlots absent. Julia 1.10.10
compatibility and all pinned manifests are unchanged; these are the actual
available Julia versions, not a claim of rerunning 1.13.2 or 1.10.10.

Coverage includes repeated/duplicate/zero factors, shared prefixes and terminals
used internally, equal operands, unused leaves, mixed signs/zero weights,
meaningful nonzero U/V, changing weights, dirty scratch, algorithm switching,
restricted-order fallback, original constructors, stale topology rejection,
pruning/empty models, custom radials, scratch growth, selected centers, generic
single-site/batch routes and serial/threaded force ownership. Existing units,
periodic multi-image/skew/open geometry, force/strain differences, owning-output
lifetimes, masks, fitting/persistence, AD force pullbacks/rollout and
feature/equivariance/parity regressions remain passing.

The new DAG force and strain finite differences retain 2e-7 tolerance at h=1e-5;
weighted-product finite differences retain the established 2e-4 Float32 / 5e-8
Float64 tolerances. Explicit reassociation observations cover arbitrary factor
orders, nonzero U/V, three input/weight changes and both references. Maximum
absolute (P, dU, dV) disagreement versus polynomial is
(1.907e-6, 4.768e-7, 4.768e-7) in Float32 and
(3.553e-15, 8.882e-16, 8.882e-16) in Float64. Against prefix it is
(2.384e-7, 2.384e-7, 2.384e-7) and
(8.882e-16, 4.441e-16, 4.441e-16), respectively. Existing 2e-5 / 3e-12 product
comparison tolerances remain unchanged. This validates reassociation error;
it does not promise bitwise equality for arbitrary factor ordering.

Direct immutable-stage-one comparisons pass **72 full E/F/V cases**: all five
fixtures at 128/1,024 atoms, selected small/large modern cases at 11,664, both
precisions and all three current algorithms. Integrated prefix is bit-identical
to immutable prefix on these outputs. DAG energy is identical there; maximum
force/virial disagreement is 2.504e-14 eV/Å / 1.819e-11 eV in Float64 and
1.336e-5 eV/Å / 0.00390625 eV in Float32. Cached independent-oracle comparisons
pass wherever references exist; all five model/oracle source identities match.
The observation table records per-model absolute errors, including larger
Float32 application sensitivities rather than concealing them behind relative
energy tolerances. No existing numerical contract was loosened.

Timing hardware/software: Ryzen 9 9950X (16 physical/32 logical CPUs), affinity
CPU 0 / 0–3 / 0–15, Julia 1.13.3, SIMD 3.8.2, SpheriCart 0.2.4,
StaticArrays 1.9.24, BenchmarkTools 1.8.0, AtomsBase 0.5.3,
AtomsCalculators 0.2.3, KernelAbstractions 0.9.45, AcceleratedKernels 0.4.4.
Exact test/benchmark manifests, hashes, LLVM identity and CPU/OS information are
retained in the new result directory. No dependency was added or upgraded by
this implementation. Real AMD hardware and optional UnicodePlots were not rerun.

## External context and limits

External comparisons use cached independent modern ACEpotentials 0.10.2 and
legacy ACEpotentials 0.6.6 / ACE1 oracles. Modern reference cache includes
128/1,024/11,664 atoms; published Ti and Ti–Al–V independent caches cover 128
atoms. Their 1,024-atom results compare all PFW algorithms, without claiming an
independent external oracle there. External Float32 execution is unsupported;
Float32 versus an independent Float64 oracle measures application sensitivity
in addition to rounded-parameter/geometry execution error. Cache hashes,
versions and timing boundaries are in the raw records. Modern public reference
constructs a fresh default PairList per call; legacy's three JuLIP calls retain
the native unchanged-geometry neighbor-list cache. Cached external timings are
context, not contemporaneous matched timings; they are not used for incremental
DAG speedup claims. No refitting, new model acquisition or PACE export occurred.

Hardware counters and cache-miss measurements are unavailable in this run. Scratch
size, static code and irregular access do not establish bandwidth saturation.
Real GPU execution was not rerun; portable KA CPU regressions are labelled as
such. This work does not change GPU, AD, feature, radial, layout or scheduling
algorithms.

## MACE interpretation

Transferable observations are deterministic structural planning, preserving
readout coordinates while sharing intermediates, additive reverse accumulation,
exclusive numerical scratch, and separately measuring construction versus
reused execution. Scalar commutative monomial sharing does not establish MACE
tensor contraction speedups, training AD speedups, or GPU performance. Whether
similar reuse outweighs indexing, scratch and preparation costs in tensor products
remains an untested hypothesis. No MACE implementation or subsequent optimization
stage is included here.

## Recommendation for review

Keep the validated DAG as an explicit opt-in for product-heavy prepared CPU ACE,
and keep `:auto` unchanged. Complete serial gains for large TiAl are repeatable
in the one-core scaling job, while small/public/parallel cases and fresh prefix
construction costs do not justify a universal promotion or a model-size heuristic
in this step. Review DAG reverse/indexing overhead as the next narrow large-model
bottleneck; treat small-model public preparation as a separate measured issue.
The current evidence does not select a cache, radial, scheduling or graph
factorization rewrite. No next-stage optimization was started.
