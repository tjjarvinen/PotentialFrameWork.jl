# GPU radial pair inference

The tested device is an AMD Radeon RX 6800 XT (`gfx1030`). GPU execution is
single-image inference only: no GPU ACE/MACE or GPU AD. Existing CPU entrypoints,
multi-image geometry and functional differentiation remain available.

KernelAbstractions and AcceleratedKernels are portable core dependencies. AMDGPU
is an application dependency: loading PotentialFrameWork does not load AMDGPU or
require an AMD driver. The real-device validation environment uses Julia 1.13.3,
KA 0.10.3, AK 0.5.10 and AMDGPU 2.10.0. AK 0.4/KA 0.9 adapters retain compatibility
with the existing SpheriCart 0.2 CPU environment. Core Julia compatibility remains
1.10.10. AMDGPU 2.8 with older KA/AK failed device compilation in the combinations
tried; this is not a supported GPU configuration. KA 0.10's POCL CPU backend and the KA 0.9 CPU backend are validated on Julia
1.10.12 and 1.13.3. An initial POCL struct-reinterpretation compilation failure
was fixed by decoding the verified species fields explicitly.
CUDA, Metal and oneAPI have not been tested.

Create an isolated application environment, from the repository directory:

```sh
julia --project=/tmp/pfw-gpu -e 'using Pkg; Pkg.develop(path=pwd()); Pkg.add(["AMDGPU", "KernelAbstractions", "AcceleratedKernels", "AtomsBase", "AtomsCalculators", "AtomsBuilder", "AtomsSystems", "Unitful", "StaticArrays"])'
julia --project=/tmp/pfw-gpu examples/gpu_pairs.jl
julia --project=/tmp/pfw-gpu test/gpu_pairs.jl --amdgpu
julia --project=/tmp/pfw-gpu benchmark/gpu_pairs.jl
```

The dedicated device test fails if the expected RX 6800 XT is unavailable; it
never substitutes a CPU run. `julia --project test/gpu_pairs.jl` runs explicitly
labelled KA CPU checks instead. Run `Pkg.test()` for CPU/interface/AD regressions;
its dependency resolution may choose KA 0.9 to accommodate SpheriCart 0.2.

## Entrypoints and ownership

```julia
using PotentialFrameWork, AMDGPU
import AtomsCalculators as AC
NL = PotentialFrameWork.NeighbourList
PP = PotentialFrameWork.PairPotentials
backend = NL.KA.get_backend(ROCArray(zeros(1)))
p = PP.lennard_jones(;epsilon=0.0103, sigma=3.4, radius=8.5)
config = NL.GPUConfig(;strategy=:central, workgroup=64, tile=4, skin=0)
w = NL.GPUWorkspace(backend; precision=Float32, config)
y = PP.evaluate_gpu!(w, sys, p; request=(AC.Energy(), AC.Forces(), AC.Virial()))
input = PP.prepare_gpu_execution(w, sys, p)
device_result = PP.execute_gpu!(input; parameters=(;epsilon=0.012, sigma=3.4))
NL.KA.synchronize(backend)
```

`evaluate_gpu(sys,p,backend;precision=Float64,config,...)` constructs a workspace;
`evaluate_gpu!` reuses it. Default computation precision is **Float64**, independent
of system precision. Float32/Float64 governs all floating point device inputs,
parameter tables, accumulation and results. Extraction explicitly converts units
to Å/eV. Host positions/system data are never modified. Every nonempty subset of
Energy/Forces/Virial, duplicates/order normalization, and empty requests are
supported. Unitful convenience results own host arrays in original atom order;
this boundary synchronizes, reads requested outputs and rejects nonfinite output.

`prepare_gpu_execution` extracts/validates geometry and uploads reusable buffers;
it performs no CPU neighbor construction. `execute_gpu!` returns owning device
arrays for requested quantities (energy is a 1×1 array), without a force readback. It reads one reduced scalar to reject singular or
nonfinite callback output explicitly, including force-only requests; this
synchronization is included in numerical timings.
It rebuilds fresh cells every call unless skin is enabled. Input upload is excluded
from this numerical call. Update `input.positions` on device in original order,
using finite coordinates in the selected precision. Reprepare after cell,
periodicity, count/order/species or cutoff changes. Preparation on the same
workspace invalidates older borrowed inputs, which subsequently reject execution.
Do not share a workspace between concurrent executions. Results remain independent
of subsequent preparations and executions. Device output copies are deliberate
ownership allocations; warmed execution is not advertised as allocation-free.

Explicit parameters refresh on every supplied call, using original model species
indices. Omitted numerical parameters retain the input's currently prepared
values. Repreparing uses the supplied model/default parameters again. SmoothPair's
triangular amplitudes and Lennard-Jones epsilon/sigma are supported. Cutoff is
structural, not a trainable execution parameter. Protocol-compatible stateless
custom radial functions can extend `PP.prepare_gpu_radial(radial,T,parameters)`;
this explicitly certifies symmetric directed semantics and selected precision.
See `DeviceQuadratic` and `SpeciesQuadratic` in the tests. General models may extend
`PP.prepare_gpu_potential` and implement `NI.evaluate_pair!` with scalar distance
and radial derivative setters. Arbitrary CPU scratch, legacy/asymmetric and
anisotropic callbacks reject GPU preparation. ChemicalSpecies accessors preserve
full named/isotopic identity; dense channels are internal parameter-table indices.
AtomsBase 0.5.3 has no integer conversion constructor for this representation:
we use its verified eight-byte isbits representation via `reinterpret(UInt64,s)`
on the host, explicitly decode its three fields on device, and test strict
`isequal` round trips, including names and isotopes.

## Geometry, strategies and skin

Only axis-aligned orthorhombic/open/mixed cells are supported. Periodic lengths
must be finite, positive and strictly greater than twice cutoff+skin. Equality,
multi-image and skew requests reject explicitly. Open axes can have zero/infinite
cell generators. Finite positions may lie outside the periodic box. Sparse sorted
integer-triple cell keys avoid vacuum-sized dense tables, with checked host bounds.
Device record sorting keeps position, species and original index aligned; device
range lookup uses binary searches over sorted occupied keys. Aliased neighboring
periodic cells are visited once. Atom identity excludes self; distinct coincident
atoms still evaluate. SmoothPair has a regular origin; singular LJ origins reject.
Callbacks require finite energies/derivatives. Invalid/singular output marks a
reduction partial; both numerical and convenience entrypoints reject it after
reading one reduced scalar.

`:central` assigns one thread per center. `:candidate` is experimental: one
workgroup owns a bounded sorted central tile; its lanes load contiguous candidate
cell runs and reuse candidate coordinates across centers in that cell tile.
Central coordinates and partial sums occupy local memory. Candidate runs larger
than a workgroup are traversed in chunks. Uniform barriers and tree reductions
include tail lanes. A center is never split across workgroups, so no inter-group
force atomics or full N-atom force buffers per workgroup are needed. Launch size
(32/64/128/256) and central tile (1/2/4) are independent of CPU SIMD width; combinations
exceeding 64 KiB local memory reject. Both strategies use the same radial callback,
full directed forces and half-weighted energy/virial. Original-order force writes
are fused into the interaction kernel. Candidate tiling may be slower.

Skin defaults to zero. Explicit `skin>0` builds exact conservative atom CSR
connectivity at cutoff+skin, after sorting. A device count/scan/fill pass sizes
connectivity without truncation. Sorting/permutation and reference coordinates
remain fixed between rebuilds; every interaction uses current positions and the
physical cutoff. Both strategies consume this connectivity. The tiled buffered
kernel follows per-center CSR lists rather than assuming stale cell adjacency.

A device displacement reduction reads one host maximum each reused evaluation.
Rebuild occurs before `2*max_displacement >= skin`, with a 16-epsilon conservative
margin. Raw-coordinate endpoint differences bound periodic torus displacement:
wrapped crossings/large jumps conservatively rebuild; winding is never inferred.
Reference positions are never advanced on a reuse step. Cell, periodicity,
count/order/species, channel order, cutoff/skin changes invalidate connectivity.
Numerical callers must reprepare metadata changes as described above. Neighbor
capacity grows on rebuild when required; numerical buffers are reused at fixed N,
and changing N reallocates them. The implementation does not tune skin.

## Measurements

`benchmark/gpu_pair_measurements.toml` records synchronized warm medians/min/max,
seven repetitions, first preparation/execution (including compilation/capacity
growth where applicable), seeds, atom counts, geometry, density, precision, model/request,
CPU/public/device workloads, separately synchronized profiling stages, host
allocation counts and device buffer sizes. The numerical workload excludes upload
and includes fresh sorting/search/reduction plus owning result copies. The public
workload includes extraction/upload and unitful host outputs. Key generation and
sort are also profiled together to sort identical unsorted records each sample;
stage timings are not additive estimates of the asynchronous total. Host output
copy includes totals and requested forces; output permutation is fused above.
Host allocation bytes are Julia wrapper allocations, not device memory allocation
counts. Device memory includes separate main, upload, algorithm scratch and owning
result buffers. No zero-allocation or cross-vendor performance claim is made.

Structures are AtomsBuilder cubic FCC argon replicated 14³ and 29³ (10,976 and
97,556 atoms), rattled within 0.03 Å using fixed seeds, then converted with
AtomsSystems.generic_system. The mixed/slab cases alternate Ar/Ne with a symmetric
SmoothPair table; slab z is open. LJ uses epsilon=0.0103 eV, sigma=3.4 Å and cutoff
8.5 Å. `gpu_skin_*.toml` records forty evolving nonuniform displacement steps,
rebuild/reuse samples, rebuild frequency and amortized time, each checked against
fresh GPU results outside timing. Compilation of both rebuild/reuse paths is
warmed before those measurements. See STATUS.md for the actual validation gates,
benchmark interpretation and remaining limitations.
