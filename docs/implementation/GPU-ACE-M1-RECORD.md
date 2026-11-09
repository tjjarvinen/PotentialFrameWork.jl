# GPU ACE milestone 1 implementation record

Implemented 2026-11-09 from 5eef1c1 under the supplied milestone instructions,
which supersede earlier GPU ACE prohibitions for this phase. The original
untracked instruction file, pinned reference environments/data, package identity
and unrelated historical artifacts are preserved. No delegation was used.

The sandbox hides /dev/kfd and /dev/dri; host-accessible rocminfo and AMDGPU
identify the RX 6800 XT, gfx1030, 72 CUs. Actual hardware validation uses Julia
1.13.3, AMDGPU 2.10.0, KA 0.10.7, AK 0.5.12 and HIP 7.1.52802. SpheriCart
0.2.4's KA 0.9 constraint requires a separate CPU oracle environment. Core Julia
compatibility remains 1.10.10, with CPU tests also run on 1.10.12. CPU-only loading
without AMDGPU passes. Device integration is an optional extension.

`adapt_model` creates independent CPU or packed device model/parameter bundles.
Precision conversion includes all numerical coefficients/configuration; masks,
indices, species and product order remain explicit. Float32 source provenance
survives promotion. Persistence keeps old checkpoints/imports and supports adapted
Float32 models. The CPU readout buffer now follows coherent coefficient precision;
functional AD entry points do not invoke inference conversion.

GPU ACE reuses sorted directed neighbor storage, physical cutoffs, skin validation
and original-index mapping. Staged density uses bounded 8×8 channel/neighbor tiles;
products/readout perform division-free reverse; the edge pass recomputes basis
values/derivatives and atomically scatters forces. U and dU remain separate, with
V/dV only when the species-order policy requires split densities. Site batches
bound scratch without losing external-neighbor forces. Results own storage by
default; borrowed execution, parameter/input updates, metadata/capacity changes
and all property subsets are tested. The local neighbor protocol remains pure.

Imported splines retain their exact piecewise cubic/transform/envelope. Custom
callables provide fixed-size values and analytic first derivatives, with explicit
configuration adaptation. Angular values/Cartesian derivatives use the existing
real-harmonic convention and are independently checked against SpheriCart.

The public Ti and April08 Ti–Al–V JSONs contain legacy precontracted complex PI
models. Acquisition retained source metadata, SHA-256 pins, CC BY 4.0 and Crown/AWE
notices. Initial sandbox DNS failure was resolved through host network access.
A pinned ACEpotentials 0.6.6 producer converts exact analytic recurrences and
algebraically expands complex products into real harmonics, retaining original
fitted coefficients/order, pair and reference terms. No refitting or spline
approximation is used. The paper's invariant counts differ from converted raw PI
coordinates; individual-row invariant inspection/selection/fitting reject these
imports. Independent original-package E/F/V checks precede timing.

Full CPU suites pass 53,958 assertions on each Julia version. The final RX ACE
suite passes 758, existing GPU pairs 1,539, and independent public CPU/AMD checks
34 per backend. Force/open-strain/periodic-strain finite differences pass the
2e-7 absolute bound at h=1e-5; detailed errors are in the validation TOML.

Quick and bounded full presets ran sequentially with synchronized completed GPU
samples, common Float64 checks, separate CPU thread processes, both precisions,
all three models, alternative mappings, fresh cells and five-frame skin/rebuild
sequences. Quick covers 11,664-atom TiAl; bounded full uses 128 atoms per model.
Larger public-model/full-default cases remain unmeasured. Compact wins strongly
at 128 atoms, wins TiAl at 1,024 and is approximately tied at 11,664. Automatic
selection is therefore limited to the RX / Float64 / small-spline 128–1,024 range;
other cases default to staged and retain overrides. The final worker records
actual/requested and complete-call mappings, rounded-CPU E/F/V errors, edges and
resident array payload categories. Historical explicit-strategy complete calls
used the previous staged default; the dated report explains that boundary.

CPU audit covers 42 TiAl and eight cases per public model, plus actual warmed
LLVM/native radial, site, reduction and SpheriCart kernels. Initial wrapper-only
radial dumps are preserved and corrected rather than interpreted as absent SIMD.
Packed Horner/compiler-vectorized SpheriCart do not imply vectorized density,
product reverse or scatter. Batched SpheriCart wins these CPU cases. Source changes
for the audit introduce no CPU optimization.

APIs/limits: `docs/gpu_ace.md`. Reproduction, attribution and schemas:
`benchmark/gpu_ace/{README,THIRD_PARTY,REPORTS}.md`. Measured CPU findings:
`docs/gpu_ace_cpu_simd_audit.md`. Raw dated reports/logs, compressed code dumps,
original oracle fixtures and rebased dependency snapshots are under
`benchmark/gpu_ace/results/2026-11-09-julia1.13.3-5eef1c1-dirty/`.

GPU training/AD, equivariant outputs, skew/multiple images and other vendors are
outside this milestone. Large legacy callable models remain expensive, especially
in the repeated-channel force pass. Device allocator peaks/register spills are
unavailable; array payloads and host allocations are reported separately. No
universal speedup, runtime SIMD percentage or Float32 suitability is claimed.
