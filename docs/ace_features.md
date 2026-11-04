# CPU equivariant ACE features

`PotentialFrameWork.ACE.ACEFeatureLayer` produces per-atom features independently
of the scalar ACE calculator. It uses the existing full directed CPU environment
and SpheriCart extension. Load `SpheriCart` to evaluate; load `ChainRulesCore` and
`Zygote` for reverse differentiation, and `ForwardDiff` for trainable callable
radial parameters or the values-only radial derivative adapter.

```julia
using PotentialFrameWork, SpheriCart
const A = PotentialFrameWork.ACE

# Radial labels specify which radial channel multiplies each input harmonic l.
# The callable includes its cutoff envelope, exactly once.
radial = A.CallableRadial(A.SmoothExponential(),
    A.RadialMetadata([(1,0), (1,1), (2,1), (1,2)]; cutoff=4.0, ad=true);
    config=(cutoff=4.0, powers=[0,0,1,0]))
layer = A.ACEFeatureLayer(
    species=[:Si, :O], radial=radial, correlation=3, lmax=2,
    cutoff=4.0, totaldegree=6,
    output_irreps=((L=0, parity=:even, channels=16),
                  (L=1, parity=:odd, channels=8),
                  (L=1, parity=:even, channels=2)))
ps = A.feature_parameters(layer;
    custom_radial=(rate=[0.4], amplitude=[1.3]))
input = A.ace_inputs(system, layer)  # system implements AtomsBase
blocks = A.ace_features(input.positions, layer, ps, input.environment, input.cell)
A.feature_basis_summary(layer)
```

The returned tuple follows `output_irreps` order. Each owning array has axes
`(magnetic component, output channel, original atom index)` with shape
`(2L+1, channels, atoms)`. Reusing inputs or a subsequent forward call cannot
overwrite a previous output. Empty and isolated systems are supported; no constant
zero-density feature is included. Species are elements, with explicit rejection of
isotope/name labels, as in scalar ACE. Numeric X is 3×N in Å, and H is 3×3 with
cell vectors in columns. `ace_inputs` converts units; open generators are zero.
Preparation is discrete and must be rebuilt when topology/species/cell changes.
The numeric function computes displacements from current X and stored image
shifts, retaining all active periodic images, including nonzero self images.

The parameter tree is
`(mixing=(block_1=W1, block_2=W2, ...), radial=(scale=...), custom_radial=...)`.
Every mixing array has axes `(coupled basis channel, output channel, central
species)`, with species order equal to the constructor. The same matrix acts on
every m component in a block. Parameters are authoritative on every call; there
is no precontracted inference cache. Default mixing repeats identity columns when
width exceeds multiplicity. `feature_basis_summary` reports multiplicity and the
maximum linear rank of each readout; additional output width does not create new
independent basis channels. Scale axes are `(radial index, center species,
neighbor species)`; callable parameters may be arbitrary nested NamedTuples of
real scalars/arrays and are shared across species. A matrix of radial callables
can supply distinct fixed central/neighbor configurations. All channel specs must
agree. Imported `RadialChannel` tables require an explicit `radial_spec` keyword
when used in the generated constructor because the tables do not store labels.

## Coupling and truncation

Density channels are species resolved: for each `(neighbor species, n, l, m)`,
sum `R_nl(r) Y_lm(r̂)` over that species' neighbors. Correlation is the maximum
number of density factors, so correlation three includes up to four-body terms
including the center. Repeated neighbor indices in products are allowed, as in
ACE density products. Input `lmax` limits factors, not intermediate or output L.
`totaldegree` bounds the sum of n+l labels in a product. Its default allows all
products through the requested correlation within the supplied radial spec.

Construction enumerates nondecreasing factor-channel lists and couples complete
multiplets along all allowed sequential angular paths. Multiplication canonicalizes
commuting density factors. Reorthogonalization within each factor list and L
removes zero/dependent paths, including self cross products. Whole multiplets
are normalized together to mean squared polynomial coefficient norm one. Parity
is `(-1)^sum(l)` and may differ from `(-1)^L`. Missing requested irreps fail,
including irreps removed by truncation. There is no feature pruning/selection API
that could cut a multiplet in half. Construction order is deterministic.

The implementation uses a construction-time factorial Clebsch–Gordan formula
with 256-bit BigFloat arithmetic and transforms to the SpheriCart real convention.
No new package dependency or core Julia compatibility change is needed.
Polynomial coefficient tables and sparse maps are Float64. Numerical evaluation
uses at least Float64 precision and preserves promoted parameter/position types
(including dual arithmetic); there is no Float32 execution mode in this phase. Reorthogonalization
rejects residual norms ≤1e-10; coefficient cleanup uses 1e-13. This is a small
CPU reference constructor, with combinatorial cost as correlation, species,
radial count and lmax grow; it is not an optimized general basis library.

Real harmonics are L2 normalized with m ordered -l:l and index l²+l+m+1.
In particular `Y_1 = sqrt(3/(4π))*(y,z,x)/r`. From complex Condon–Shortley
harmonics C, for m>0,
`Y_(+m) = ((-1)^m C_(+m) + C_(-m))/sqrt(2)` and
`Y_(-m) = ((-1)^m C_(+m) - C_(-m))/(i sqrt(2))`.
The real coupling includes phase `i^(l1+l2-L)`; the distinct-vector `1o × 1o → 1e`
channel is `-cross(a,b)/sqrt(2)` in the cyclic `(y,z,x)` ordering. Rotations act
on column m components through the real harmonic D_L matrix, for both parities;
inversion multiplies each block by its declared parity. The tests compare these
statements to direct Cartesian cross products and independently sampled harmonic
rotation matrices through L=3.

## Resolved reference tables and scalar compatibility

For hand-resolved multiplets use
`ACEFeatureLayer(radial, density, products, coupling; species, radial_spec, lmax,
output_irreps, cutoff=nothing)`. Density entries are `(neighbor channel, radial
index, harmonic index)`; neighbor channel 0 sums species. Products hold density
indices. Each sparse coupling matrix has `(2L+1)*basis_channels` rows with m
fastest, and one column per product. The constructor checks dimensions/indices;
the caller supplies valid equivariant intertwiners in the stated convention.
The generated constructor handles coupling validity and symmetry automatically.

`ACEFeatureLayer(resolved_model; channels=1)` keeps the identical resolved scalar
0e coordinates, including row masks. Set `ps.mixing.block_1[:,1,:] .= scalar_ps.WB`
to reproduce the scalar **many-body term** summed over sites. Reference energies
and separate pair terms remain in `ace_evaluate` and must be accounted for
separately. This adapter rejects species correlation-order projected models;
it preserves the existing scalar import/native checkpoint formats and APIs.
It does not introduce a new feature checkpoint format or scalar fitting route.

## Differentiation and Lux

`ace_features(X, layer, ps, environment, H)` has a first-order ChainRules rule
for positions and explicit parameters, tested through Zygote with arbitrary output
cotangents, zero blocks and thunked blocks. The reverse contracts mixing gradients
over m, applies transposed mixing/coupling maps, and uses the existing division-free
product reverse followed by radial/angular reverse and endpoint scatter. It never
forms a full position/parameter Jacobian. Pullbacks own snapshots of numerical
inputs, model tables, parameters and topology; they recompute local forward data
from those snapshots and survive later forward calls and caller mutation.
Cell derivatives explicitly report unsupported, and topology/model data are
nondifferentiable. Fixed topology must remain valid in a neighborhood of X.

Custom radial callables use the existing
`(values, derivatives_or_nothing, r, ps, config)` contract. Analytic d/dr must match
values. Parameter derivatives require generic values code compatible with
ForwardDiff and an `ad=true` declaration; the fallback differentiates only the
local contracted radial values, not the whole feature calculation per parameter.
The existing `derivative=:ad` values-only adapter is also supported. Fixed cutoff
and envelopes are not trainable. Zero density factors are valid; coincident
physical atoms remain unsupported by this harmonic path.

```julia
using Zygote, LuxCore
energy(x, p) = sum(sum(abs2, block) for block in
    A.ace_features(x, layer, p, input.environment, input.cell))
forces = -Zygote.gradient(x -> energy(x, ps), input.positions)[1]

lux = lux_layer(layer; parameters=ps)
parameters, state = LuxCore.setup(LuxCore.Random.Xoshiro(7), lux)
features, newstate = LuxCore.apply(lux, input, parameters, state)
```

The Lux wrapper follows the existing scalar/pair adapters, retains explicit
parameters and passes state through without putting scratch in state. The
invariant head above is a composition example, not a new site-energy API.
Only first derivatives through the feature boundary are claimed. Force-loss
training and mixed derivatives through composed models need a separate validation
milestone. Feature SIMD optimization, GPU execution, full MACE interactions and
general equivariant training remain outside this phase. Scalar ACEfit is unchanged.

## Reproduction

In an application environment with this checkout developed and the test extras
installed, run `julia --project=<env> -t4 test/ace_features_runtests.jl` for the
focused numerical/derivative checks, or `test/runtests.jl` for all regressions.
`julia --project=<env> -t4 benchmark/ace_features.jl` reports warmed serial
allocations, output/saved-pullback storage and small-case timings. Topology
preparation is outside these numerical measurements. These reference paths
allocate; timings are diagnostics, not a performance claim or a regression gate.
