# Scalar GPU ACE inference

GPU ACE is explicit inference. It does not provide GPU AD, GPU features or MACE.
Existing CPU functional AD entry points retain their behavior. AMDGPU is optional;
loading the package does not load a vendor runtime. Use an application environment
with AMDGPU 2.10 / KernelAbstractions 0.10 for the tested Julia 1.13 device stack.
SpheriCart 0.2.4 constrains KA to 0.9, so independent SpheriCart checks run in a
separate CPU environment. Core compatibility remains Julia 1.10.10.

```julia
using PotentialFrameWork, AMDGPU
const A = PotentialFrameWork.ACE
model, ps, state = A.load_resolved("data/ace_ref/resolved_model_v1.toml")
bundle = A.adapt_model(model, ps; array_type=ROCArray, precision=Float64)
executor = A.GPUExecutor(system, bundle.model, bundle.parameters;
    strategy=:staged, batch_size=256)
y = A.execute_gpu!(executor)       # owning device energy/forces/virial
host = A.gpu_result(y)             # numeric Å/eV arrays, original atom order
A.prepare_executor!(executor, next_system) # explicit upload; cell/species invalidate topology
A.execute_gpu!(executor; parameters=new_device_parameters)
```

`load_resolved` and `load_recipe` accept the same `array_type` and `precision`
keywords. `adapt_model` returns independent model and parameter storage. All
floating model tables and explicit floating parameter arrays use the requested
precision; integer indices stay integer. Float32 → Float64 cannot recover lost
coefficient precision. Custom config types can extend `A.convert_data(config,T,A)`;
ordinary nested NamedTuples/tuples and arrays are handled recursively. Avoid using
conversion inside differentiable functional CPU calls.

`ACECalculator(model, ps; array_type=ROCArray, precision=Float32)` supplies a
convenience AtomsCalculators wrapper. Its complete calls include extraction,
upload, preparation, execution and owning unitful host results. The numerical
executor accepts `properties=(AtomsCalculators.Energy(), ...)`, including all
seven nonempty subsets and `()` (returns an empty NamedTuple and resets outputs).
`owning=false` returns borrowed numerical output views, invalidated by the next
execution or preparation. Ordinary results remain valid after later calls/growth.
Scalar and resident array parameters, including nested radial scales and custom
radial parameters, support explicit updates; scalar types and array shapes must
remain fixed. Type/shape changes require adaptation and preparation.

GPU geometry reuses NeighbourList sorted sparse cells and complete directed CSR.
Periodic lengths must strictly exceed twice cutoff + skin. Orthorhombic, open and
partially periodic cells are supported; skew and multiple-image GPU geometry are
rejected. CPU multi-image/skew support is unchanged. Central selections use original
atom indices; neighbor forces may reach atoms outside the selection/site batch.
The physical radial cutoffs remain active with skin candidates. Coincident physical
atoms reject, matching existing scalar ACE behavior. Direct mutation of resident
positions requires finite coordinates and unchanged cell/species; use preparation
for metadata or capacity changes. Positive skin uses the existing conservative
maximum-displacement validity check; zero skin rebuilds cells/connectivity each call.

Staged execution uses a 64-thread workgroup per site with bounded eight-channel
by eight-neighbor tiles and cooperative density reductions, completes
densities at a kernel boundary, contracts products/readout with division-free
polynomial reverse, and revisits edges for analytic forces/virial. Separate U/V
and adjoint arrays retain species-order projection; V and its adjoint have zero
size when the selected order range needs no split density. No per-edge basis tables or
full-system per-worker force arrays are stored. The central force is locally
summed, and both endpoints use atomic addition. W = -Σ displacement ⊗ gradient.
Parameter precontraction refreshes every execution from explicit readout arrays.
`:baseline` uses one site/thread in the density pass. Experimental `:compact`
uses four logical 16-lane site groups per workgroup. Each site owner contracts
products; a barrier publishes adjoints before cooperative edge revisits and local
energy/central-force/virial reduction. Density remains a separate cooperative pass.
Only bounded reduction tiles occupy shared memory; no hardware warp width is assumed. `:auto` selects compact at initial executor preparation only on the measured RX
6800 XT / Float64 / 128–1,024-atom range for spline-only models with lmax ≤ 2,
at most 27 density components and 167 products. It selects staged otherwise.
These conservative dimension/size bounds come from TiAl measurements; they do
not assert a universal speedup. Float32 and large/custom models retain manual
compact selection pending broader calibration. Construct a new executor to
reselect the strategy for a different workload; per-step switching is absent. `gpu_memory(executor)` reports
array payloads separately from host allocations and auxiliary backend workspaces.

Imported cubics retain exact spline tables, transform, envelope and first derivative.
Custom analytic callables use the same in-place radial protocol with fixed output
count; GPU execution uses fixed-size temporary vectors. Values-only AD wrappers
and allocating/host-only operations are not supported on device. Compilation errors
identify unsupported device operations; CPU-only callables remain usable on CPU.
User code owns smoothness, branch costs and valid device operations. Custom
configuration structs need both precision conversion and device argument
adaptation (for example `Adapt.@adapt_structure`); NamedTuples work recursively. Registered
legacy orthogonal polynomials are analytic; they are not spline approximations.

`compare_precision(model, ps, system_or_collection; array_type=Array,
properties=..., tolerances=(forces=...))` returns per-structure and pooled reports.
No tolerances gives descriptive status. Tolerances check absolute total energy
and maximum force/virial component error. RMS pools all components, not structure
RMS averages. Failures/nonfinite values stay visible and make aggregate errors
infinite. Empty collections report `:empty`. Units, source precision, geometry
rounding, component/atom identifiers and model/structure hashes are retained.
This measures supplied examples, not universal Float32 suitability.

Public Ti/Ti–Al–V acquisition, original-package oracle and benchmark commands are
in `benchmark/gpu_ace/README.md`. Legacy JSONs contain precontracted complex PI
coordinates. Conversion preserves those explicit coefficients and expands their
real parts algebraically into real harmonics. These coordinates are not the
paper's 1,809 / 29,703 independent invariant fitting coordinates. Individual-row
inspection/selection and invariant fitting reject this imported representation;
full scalar inference, central selection and species-order projection remain
mathematically defined. The original source files, hashes and ordering are retained.

Measured timings and limits are in the dated benchmark results; CPU source/native
evidence is in [the SIMD audit](gpu_ace_cpu_simd_audit.md).
