# CPU ACE validation and measurements

The reference environment stays in `data/ace_ref`, pinned to Julia 1.13.2. The
core package retains Julia compat 1.10.10 and is tested on 1.10.12 and 1.13.2.
Do not run `generate_reference_data.jl` to validate: it fits and overwrites the
saved model. These scripts only load the saved model:

```sh
julia --project=data/ace_ref data/ace_ref/validate_saved_model.jl
julia --project=data/ace_ref data/ace_ref/export_resolved_model.jl
julia --project=data/ace_ref data/ace_ref/validate_resolved_model.jl
julia --project=data/ace_ref data/ace_ref/validate_derivative_products.jl
julia --project=data/ace_ref --threads=4 data/ace_ref/validate_cpu_model.jl
julia --project -e 'using Pkg; Pkg.test(;julia_args=["--threads=4"])'
```

The supplied JSON is a recipe, not an explicit basis dump. The producer exports
the constructed spline radial/pair tables, Agnesi transforms, envelopes, density
indices, selected products, sparse invariant map and original parameter ordering
to `resolved_model_v1.toml`. Its 270 parameters are 246 ACE readout values and 24
pair readout values. Spline coefficients use a numerical change of cubic
coordinates, checked at independent interior points. Replacing these splines
with their generating polynomial recurrence changes the reconstructed model.

`reloaded_predictions.toml` records full precision outputs on the rounded XYZ
geometries. Tests compare residual/component-sensitive absolute tolerances, not
only relative total energy dominated by one-body references. The strict element
channel policy rejects named species and isotopes. Nonzero self-images and
multiple images are retained; only the zero-shift self edge is excluded.

The serial implementation materializes radial/angular/density/product data for
diagnostics. The prepared executor capacity-manages numerical buffers. Every
centre's advertised footprint contains the centre and every full directed
neighbour endpoint. Phase tasks use exclusive numbered scratch slots independent
of thread IDs and accumulate into one shared species-sorted 3×N force array.
Owning outputs, saved pullback copies and local edge/feature buffers are separate
memory costs. Preparation still allocates snapshots and topology; there is no
Verlet reuse. `prepare_executor!` follows every environment rebuild, including
same-sized changed species/graphs. Precontraction refreshes from current ps.

Linear features preserve original `vec(WB),vec(Wpair)` coordinates. Rows are E,
`vec(F)`, `vec(W)`; baseline subtraction is explicit. `training_rows` omits missing
properties/`nothing` labels. Weights multiply residuals directly. ACEfit 0.3.1
uses `solve(QR(lambda=..., P=prior), A, y)["C"]`; the prior enters as lambda*P in
the augmented solve. Smoothness is the documented metadata-derived algebraic
prior in `smoothness_prior`, applied to original coefficients. It is not claimed
to reproduce every ACEpotentials prior. `fit_pair=false` holds the imported pair
readout fixed; default fitting solves ACE/pair columns jointly. Rank-deficient
tests recover predictions rather than asserting unique coefficients.

Nonlinear radial/pair scales are explicit nested parameters from
`trainable_parameters`. The resolved spline tables/cutoffs stay fixed. LuxCore
layers expose the same parameter tree and fixed state; gradients through force
outputs are implemented at the numerical primitive. Position force pullbacks use
a directional derivative of the analytic reverse (Hessian action); mixed products
reuse that direction in analytic parameter gradients. No full Hessian/Jacobian or
whole-system replay per parameter is used. Pullbacks deep-copy inputs/maps/model
data for a safe lifetime. Virial position products are supported; cell derivatives
explicitly report unsupported requests. First reverse derivatives through forces,
Zygote, local ForwardDiff references and short fixed-cell rollouts are tested.
Arbitrary higher reverse orders and AD through mutable inference/preparation are
outside the verified interface. See [the training example](../../examples/ace_training.jl).

Run measurements without other benchmark/test jobs:

```sh
julia --project=data/ace_ref --threads=1 benchmark/ace/run.jl
julia --project=data/ace_ref --threads=4 benchmark/ace/run.jl
julia --project=data/ace_ref benchmark/ace/inspect_simd.jl
```

Reports include all samples, extrema/median, allocation bytes, BLAS/Julia threads,
versions, CPU, preparation stages, phase occupancy, species packing and numerical
scratch/force storage. They measure two supplied fixtures and synthetic 8/64-atom
periodic systems with one/two populated element channels and two spacings. BLAS
uses one thread; affinity is not pinned. Oracle timings use the same saved model
and geometry with its public default calculator. Both public E/F/W timings include
list preparation and unit/output conversion. Prepared timings exclude preparation
and owning outputs and are labelled separately. AD save/backward measurements
cover the two-atom periodic fixture only. First-site stage timings are diagnostics,
not an additive decomposition of a full parallel run.

Recorded on AMD Ryzen 9 9950X, Julia 1.13.2, BLAS one thread, width 4, default
SpheriCart angular backend (nine samples after three warmups). Public E/F/W
median times include preparation, units and owning outputs, in microseconds:

| Fixture | Atoms | Julia threads | Independent CPU | ACEpotentials 0.10.2 |
| --- | ---: | ---: | ---: | ---: |
| 1 | 2 | 1 | 51.920 | 81.150 |
| 1 | 2 | 4 | 52.060 | 84.390 |
| 33 | 54 | 1 | 1144.383 | 1599.675 |
| 33 | 54 | 4 | 1151.104 | 736.322 |

These cases demonstrate the importance of comparing thread settings and complete
calls. The full-site write footprints limit available parallelism in dense
fixtures. The four-thread prepared baseline allocates 64 bytes for fixture 1
and 1,728 bytes for fixture 33; owning public calls allocate 79,144 / 662,392
bytes respectively. The shared force array occupies 48 / 1,296 bytes at both
thread counts. Site scratch occupies 57,448 bytes with one slot and 229,672 bytes
with four slots on these fixtures; it stores local numerical features, not four
full-system force arrays. Preparation remains a substantial allocation source.

The final two-atom force-loss pullback takes 126.691 microseconds median with
625,424 allocated bytes in the four-thread report; its owned saved state occupies
1,135,889 bytes. These are measured correctness-oriented AD costs, not a training
throughput guarantee. The superseded prototype used a repeated local-gradient
position product and is retained only as labelled historical evidence.

Angular measurements cover lmax 1/2/4/6, 16/64/256 neighbours, widths 2/4/8 and
value/value-gradient workloads. Unit checks also cover lmax=0 and widths 1/16.
SpheriCart 0.2.4's generated and batched gradient kernels assume lmax>=1; the
adapter supplies the exact constant/zero derivative at lmax=0. The explicit
SIMD block recurrence follows the inspected L2 normalization/order and carries
Cartesian first derivatives through normalization. The block value-only path
currently computes joint values/gradients internally. Its recurrence temporaries
are invocation-local allocations. Tiled radial-angular-density fusion retains
intermediates for the reverse force contraction; it does not claim register
residency or a fully allocation-free fused kernel. Default SpheriCart remains
selectable and selected; generated/block/fused variants are explicit experiments.

For lmax=2, 64 neighbours, width 4, value/gradient output and four Julia threads,
standalone medians are 0.700 microseconds for SpheriCart batched, 1.530 for its
generated single-input loop and 3.340 for explicit blocks (48 / 0 / 2,256 bytes).
The generated loop uses a concrete-type function barrier. Complete force timings
and the other cutoffs/counts are in the reports; standalone rankings do not select
an inference backend automatically. At width 4 with SpheriCart and four threads,
tiled fusion is slower on both supplied fixtures in this report and stays optional.

Float32 executors store Float32 buffers with mixed Float64 imported radial
coefficients. Radial transform/envelope derivatives and spline interval selection
are scalar per lane; the joint cubic Horner arithmetic is SIMD. Thread/task
creation and some angular scratch allocate; public
preparation and AD saved state allocate as reported. No universal speedup,
automatic crossover, cell/pressure-coupled MD derivative or GPU support is claimed.
SpheriCart itself declares KernelAbstractions/GPUArraysCore; this package adds no
GPU execution backend, CUDA/AMDGPU dependency or device implementation.

Formula/convention adaptations and derived tables are attributed to ACEpotentials
0.10.2 and SpheriCart 0.2.4; their MIT notices are retained under
[`docs/implementation/licenses`](../../docs/implementation/licenses).
EquivariantTensors 0.4.3 is used by the producer; runtime products/couplings use
exported numerical indices and standard SparseArrays. Historical first-phase pair
benchmarks are unchanged. `results-prototype-position-vjp-*.toml` is superseded
prototype evidence, not the final directional-Hessian implementation.
