# Serial numerical site derivatives

`NeighbourInterface` (NI) declares the local protocol. `NeighbourList` (NL)
implements parameter buffers, CPU batch fallbacks, chunk reduction and the
numerical AD boundary. Models implement only local mathematics.

```julia
outputs = NI.SiteOutputs(energy=true, gradients=true, features=false)
y = NL.site_evaluate(X, prepared, ps, neighbors; outputs, centers=nothing,
                    scatter_capacity=4096, batch_size=1)
# X: 3×N numeric Å, ps: explicit nested parameters, neighbors: prepare_neighbors
# y: (energy, forces::3×N, virial::3×3, features)
```

The features field is `nothing` when absent, otherwise an original-order vector
of model-defined tuples of numeric arrays. Local blocks have no prescribed
magnetic/channel axes. ACE exports `(m,channel,original atom)` arrays through
`ace_features_chunks`; `ace_evaluate_chunks` exports the usual scalar E/F/W tuple.
Prepared neighbors must use the model's species ordering. Numeric calls recompute
all active image displacements from X and fixed shifts/cell; stored forward
positions are never substituted for X. Topology must remain valid nearby.

A model author implements:

- `NI.make_scratch(prepared, max_padded_edges, T)` for exclusive site scratch.
- `NI.make_site_output(prepared, max_padded_edges, T, outputs)` and
  `output_energy`, `output_gradient(out,k)`, `output_features`. The latter returns
  the model-defined tuple of numeric block arrays. Unrequested accessors need
  not be mathematically evaluated. Empty tuples cover models without features.
- `NI.evaluate_site!(out, prepared, site, ps, scratch, outputs::SiteOutputs)`:
  overwrite requested outputs. `g[k]` is ∂E/∂r_k. Skip inactive edges before any
  invalid geometry operation. Coupling between neighbors is allowed.
- `NI.site_derivative_capabilities(prepared)` returning
  `(energy=SiteInputs(...), gradients=SiteInputs(...), features=SiteInputs(...))`.
  Each `SiteInputs(displacements,parameters)` separately declares support.
- `NI.site_vjp!(rbar,pbar,prepared,site,ps,bars,scratch; request=SiteInputs())`:
  **add** local derivatives to the requested destinations and return `nothing`.
  `rbar` is a one-based padded-edge vector of numerical three-vectors;
  `bars::SiteCotangents` has `energy`, `gradients` and `features`. Each field may
  be `nothing` for mathematical zero. `gradients[k]` uses the same padded edge
  axis; features retain the model-defined local block axes. All arguments and
  buffers must be nonaliasing. Borrowed views and scratch cannot escape.

Parameter buffers allocate once per pass, reset to zero, and accumulate serially.
Default `NI.parameter_buffer(prepared, ps, T)` preserves nested NamedTuple names,
array shapes and scalar leaves, storing scalars as mutable `Ref`s. Add a standard
parameter tree with `NI.add_parameter_buffer!(prepared,buffer,tree)`; obtain an
owning standard tree with `NI.parameter_value(prepared,buffer)`.
`NI.reset_parameter_buffer!` resets every leaf. Empty trees are valid. Models
can specialize these hooks without exposing a flat vector. Callbacks never mutate
ps. No per-site parameter gradients are retained as a tape.

For batches, `NI.evaluate_sites!` and `NI.sites_vjp!` have serial active-site loop
fallbacks using `site_count`, `site_at` and `output_at`. Numerical forward batches
are bounded by `batch_size`; the backward executor currently calls `site_vjp!`
directly. Batched inference alone does not establish batched differentiation.

`NL.site_pullback(X,prepared,ps,neighbors,bars; request=SiteInputs(), ...)` exposes
the numerical backward pass for independent checks and parameter-only requests.
Global bars are `(energy,forces,virial,features)`; absent fields are `nothing`.
For each active edge it supplies
`gbar = Fbar_i-Fbar_j-transpose(Wbar)*r` to the model, then adds `-Wbar*g` to the
returned displacement bar and scatters minus at i, plus at j. Frobenius matrix
contraction is used. Repeated images retain separate geometry; self-image
endpoint position contributions cancel. Selected centers may affect unselected
neighbors. Force/feature cotangents are permuted exactly once to cell order.
Unsupported active combinations fail before model callbacks. Zero cotangents
never activate an unsupported capability.

The ChainRules rule for `site_evaluate` normalizes zeros/thunks before model
callbacks and snapshots X, parameters, prepared model tables and neighbor
metadata, including fixed cell and maps. It saves no workgroup or worker scratch.
Backward recomputes local intermediates, supports delayed/repeated calls and
returns owning gradients. Saved input/topology storage scales with the system;
reusable scatter storage is bounded by chunk capacity. The serial backward
uses one global parameter accumulator and one final position array. Allocating
local reference mathematics are distinct from bounded scatter storage.

The numerical boundary is first reverse order, including local second/mixed
products for models explicitly declaring gradient-output derivatives. No arbitrary
nested reverse AD, trainable topology/cell/species, incoming-feature derivatives,
threaded backward or GPU AD is implied. Feature precision remains Float64 or
promoted arithmetic; scalar precision follows explicit model/input arithmetic.

For Lux, `lux_layer(model; parameters=ps, executor=:chunks)` consumes an input
with `positions`, `neighbors`, optional `centers`, `batch_size` and
`scatter_capacity`. State is passed through unchanged. The existing default
functional Lux routes remain available. Numerical parameters remain authoritative
on every call; derived inference coefficients do not enter this derivative path.

`input_derivatives=NI.SiteInputs(displacements=false, parameters=true)` on the
numerical forward/rrule boundary explicitly suppresses unneeded input derivatives;
the unitful explicit-parameter calculator rule uses this mode. Ordinary reverse
AD cannot infer input activity from a callback's captures, so numeric/Lux callers
can select it when training only parameters.

The optional four-argument allocation hook `NI.make_scratch(prepared,capacity,T,ps)`
is used by backward execution; its default calls the existing three-argument hook.
ACE uses it for one reusable ordinary and one directional parameter-reverse tree
per exclusive pass. Custom radial parameter contraction reduces the complete
site's radial contributions before local parameter AD. No full parameter tree is
allocated/retained for each site in the production chunk adapters. Local geometry,
angular/radial/product and small head-AD temporaries still allocate; these reference
mathematics are not claimed allocation-free. `NI.site_numeric_type(prepared,T)`
(default identity) can preserve a model's minimum arithmetic precision.

For a single borrowed site (typically `NI.site_at(batch,1)`), a local check is:

```julia
request = NI.SiteInputs(displacements=true, parameters=true)
capacity = NI.neighbour_count(site)
scratch = NI.make_scratch(prepared,capacity,Float64,ps)
rbar = zeros(SVector{3,Float64},capacity)
pbar = NI.parameter_buffer(prepared,ps,Float64)
bars = NI.SiteCotangents(1.0,nothing,nothing) # energy-only VJP
NI.site_vjp!(rbar,pbar,prepared,site,ps,bars,scratch; request)
parameter_gradient = NI.parameter_value(prepared,pbar)
```

A batch uses the same accumulator and one exclusive scratch; each local
three-vector destination is independent and initialized by the caller:

```julia
NI.reset_parameter_buffer!(prepared,pbar)
rbatch = [zeros(SVector{3,Float64},NI.neighbour_count(NI.site_at(batch,k)))
          for k=1:NI.site_count(batch)]
batch_bars = [NI.SiteCotangents(1.0,nothing,nothing)
              for k=1:NI.site_count(batch)]
NI.sites_vjp!(rbatch,pbar,prepared,batch,ps,batch_bars,scratch; request)
```

The fallback visits only `site_count(batch)` active sites, including a partial
last batch. Sized buffers must cover each site's padded edges. In production,
keep batch size bounded and reuse these arrays. The independent model in
[test/site_pullbacks_runtests.jl](../test/site_pullbacks_runtests.jl) supplies
complete analytic callbacks for `q=sum(abs2,r_j)` and
`E=a^2*q/2+b*q^2/4`, including a feature map `(q,a*q^2)`. Its nonlinear cross-neighbor
term prevents an edge-independent derivative implementation from passing.

Explicit-parameter `AtomsCalculators.calculate(properties,sys,calc,ps[,state])`
uses the same numerical rule for scalar ACE parameter training, with a fresh
mathematical geometry snapshot. It preserves state identity and its active-call
and failure guard; it does not save or borrow that state's inference scratch.
Input-system differentiation stays at the explicit numerical X boundary. Existing
inference state methods, colored execution, checkpoint formats and default
functional numerical/Lux paths remain available.
