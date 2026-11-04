# CPU equivariant ACE feature layer

## Starting point

Implement in PotentialFrameWork.jl from main. The inspected main revision is `13f5cf5e249692bbe0cdb0890d1476b09cdaf339`. Check the actual checkout and read the repository's AGENTS.md before editing. Preserve local changes and reference data. This task extends CPU ACE with a per-atom feature layer that can serve as the first interaction stage of a future MACE implementation.

Start with `src/submodules/ACE/{ACE,resolved,serial,projection,radials,cpu}.jl`, the ChainRulesCore, ForwardDiff, LuxCore and SpheriCart extensions, and `test/ace_cpu_ad.jl`. In the existing serial evaluator, site_data forms density and product terms, then computes B and a scalar energy. site_reverse seeds the product reverse with weights from the scalar readout. Use these boundaries to identify reusable code. The prepared CPU executor precontracts scalar readout weights into the product basis.

## 1. Define the interface

Use a separate ACEFeatureLayer entry point. Follow existing naming and optional-extension conventions, and document the concrete construction and evaluation signatures. The intended constructor contract is equivalent to

```julia
ACEFeatureLayer(
    species = [:Si, :O],
    correlation = 3,
    lmax = 2,
    # Supply the radial basis, cutoff and basis truncation settings.
    output_irreps = (
        (L=0, parity=:even, channels=16),
        (L=1, parity=:odd,  channels=8),
        (L=1, parity=:even, channels=2),
    ),
)
```

- Correlation is the maximum number of density factors. A correlation of three allows terms of up to four-body order including the center. Input lmax limits the spherical harmonics used in each density factor. Requested output L is a separate quantity governed by the allowed couplings and the basis truncation.
- Accept the explicit named tuples. Validate integer L ≥ 0, positive integer channel counts, supported parity symbols, and unique (L, parity) pairs. Keep caller order and retain typed irrep metadata alongside the layer.
- Return a tuple of owning numerical arrays in that order. Each shape is (2L + 1, channels, atoms), with magnetic components ordered -L through L in the documented real harmonic convention and atoms in original input order. Use the existing numeric position and cell conventions at the functional boundary.
- Return per-atom features. Keep the scalar energy head separate. Any future GPU execution may use sorted atom order internally with an explicit shared permutation between layers.
- Keep model and basis data, explicit parameters, topology, and scratch separate. Provide a LuxCore setup/apply route consistent with the current extensions. Keep readout weights trainable and use the same weight matrix for every magnetic component within each block. If weights depend on central species, use an explicit parameter axis and document it.

## 2. Establish the serial reference

First implement and test the feature evaluation and pullback using small resolved coupling tables. Include 0e, 1o, and a correlation-2 1e block obtained by coupling two distinct 1o density channels. A self-cross-product vanishes, so use distinct radial or species channels and a geometry that gives a nonzero result. Check this fixture against a direct cross product after accounting for harmonic convention and normalization.

Reuse radial and angular evaluation, product evaluation, and the existing full directed environment where practical. Extract shared density and product stages only where it reduces duplication without changing the scalar path. Keep species-dependent ACE density construction for this phase. Make coupling consume density/product data independently of how species or future learned embeddings produced the channels.

Implement the site calculation as products → fixed coupling into complete irrep multiplets → trainable channel mixing. Keep a simple numerical reference path. Preserve existing scalar calculator behavior, reference and pair terms, imported and native model formats, prepared CPU optimizations, and the scalar ACEfit route.

## 3. Implement the feature pullback

Add a ChainRules pullback for the functional CPU feature entry point that accepts arbitrary cotangents in the output layout, including zero and thunked blocks. For each site, accumulate the mixing-weight gradient over magnetic components, apply the transposed mixing map and coupling map to form product adjoints, then reuse product and density reverse logic. Accumulate position gradients at the center and neighbors and return parameter gradients in their original structure. Include radial parameters supported by the feature path. Document the derivative requirements for any supported custom radial callable.

Compute displacements from the supplied positions and stored periodic shifts at fixed cell and topology. Follow existing conventions for unsupported cell or topology derivatives. Preserve the forward values needed by the pullback across scratch reuse and subsequent parameter updates. Avoid materializing full position or parameter Jacobians. Check a composed invariant scalar energy head with forces obtained through the feature pullback.

## 4. Build the basis from user settings

The resolved-table prototype is an intermediate milestone. Complete the user-facing constructor by generating the density/product specification and sparse coupling for the requested irreps from the supplied radial, species, correlation, lmax, cutoff and truncation settings. Reuse suitable existing coupling machinery where available. Verify the real-harmonic conventions and record any dependency and compatibility implications.

- Enumerate the allowed products under the basis truncation. Use angular coupling selection rules and product parity, and account for symmetry of commuting density factors. Generate whole multiplets and handle zero or dependent coupling paths. Keep ordering deterministic.
- Distinguish basis multiplicity from requested readout width. Report the number of available coupled basis channels by irrep. Reject a requested block when no basis channels survive. Allow a readout wider than the basis and document its linear rank limit.
- Validate generated tables as well as the hand-checked fixture. Exercise correlation three and a higher-L output, in addition to a case with both parities at the same L. Ensure any pruning or selection exposed by the feature API acts on whole multiplets.

## Acceptance checks

- Confirm block shapes, ordering, original atom order and validation of invalid or unavailable irreps. Test rotations, inversion, translation and atom permutation using a nonsymmetric geometry. Use nonzero outputs and compare rotations in the documented real harmonic convention.
- Reproduce the existing scalar ACE many-body term using the identical resolved 0e basis and weights. Compare against the scalar path with the reference and separate pair contributions accounted for explicitly.
- Compare position and parameter directional derivatives against central finite differences for the contraction of the output with arbitrary fixed cotangents. Include a non-scalar block and a single active block, a multispecies case, and periodic geometry with fixed valid topology. Cover an isolated atom and check pullback lifetime after a later forward call.
- Run the relevant existing CPU ACE, import, persistence, training and AD regressions. Update the implementation status with commands, actual results and any remaining limitations. Use the serial CPU results as the correctness reference for later GPU planning.

## Scope

This phase covers the CPU feature contract, basis construction, and first-order differentiation into positions and explicit parameters. Schedule feature-path SIMD optimization, GPU kernels, full MACE integration and general equivariant training separately. End-to-end force-loss training needs mixed position–parameter derivatives through the composed model and a dedicated validation milestone. The existing scalar ACEfit workflow remains the fitting interface for scalar ACE.
