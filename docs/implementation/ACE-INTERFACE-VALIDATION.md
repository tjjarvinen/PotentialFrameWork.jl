# Validation requirements

Use independent mathematical checks, not tests that only replay implementation
logic. Keep existing CPU ACE, pair, image geometry, AD and persistence tests. Run
focused tests during development and full suites after integration. No brittle
timing assertions. Historical results in STATUS are context, not new evidence.

## Inspection, masking and pruning

- Reconstruct metadata from real coupling/product dependencies and compare against
  imported specs; include a row combining multiple products/orders.
- IDs/provenance remain traceable through two successive reductions and checkpoint
  round trips. Reject wrong-basis selectors and invalid indices clearly.
- Compare masked full basis versus structurally pruned E/F/W, with random nonzero
  parameters, nonunit nonlinear scales and parameter updates. Include mask removal,
  all/none selection, shared dependencies and zero density factors.
- Optimizers/readout installation cannot reactivate disabled entries. Masked
  parameter sensitivities are zero in original coordinates. Verify coefficient
  transfer and pullbacks through mappings; input models/parameters are not mutated.
- Reduced structural dimensions/dependency maps are consistent and unused work is
  removed where claimed. Changing learned zeros must not change structural identity.
- Restricted construction reaches the same retained evaluator without instantiating
  the full executor. Test old import defaults and unsupported metadata criteria.

## Center subsets and force accounting

- Default/all/range/noncontiguous/reordered/empty selections; invalid and duplicate
  indices; original versus internally sorted positions. A calculator reused after
  changes does not retain a stale center list or stale results.
- Selected E equals an independent sum of existing site_diagnostics energies.
  Subset partition sums reconstruct full E/F/W under the site's actual convention.
- At least one fixture has a selected center whose nonselected neighbor receives
  a nonzero force. Compare forces with finite differences of selected energy into
  ALL positions. Check net force, translation/rotation and strain virial convention.
- Include tiny skew periodic cells, repeated atom indices at different image shifts,
  nonzero self-images, empty species groups and all CPU SIMD tails.
- Selected phase tasks retain disjoint full force destinations; compare scalar,
  serial CPU and threaded paths. Full-system force memory must not grow N*workers.

## Species-order algebra

- For tiny products, independently enumerate factor choices in U+tV and compare
  EACH k. Include repeated factors, zeros and products of differing orders.
- Sum k=0:nu reconstructs original B and E/F/W; k=0 equals evaluation excluding
  selected neighbors with the same model, not refitted radial functions.
- Sum k>=1 equals full minus excluded E/F/W. A simple quadratic fixture shows
  mixed selected/rest terms survive k=1 while selected-selected terms do not.
- Center species must not alter selected-neighbor factor counts. Multiple selected
  species share one count; order means factors, not distinct atoms.
- No selected neighbors, all selected neighbors, absent known species, malformed
  ranges, maximum above total order, and defined empty-group behavior.
- Reference contributes at k=0; pair edges contribute at k=0 or 1 by neighbor
  species. Test with nonzero reference/pair coefficients so accidental cancellation
  or omission cannot pass. Document directed pair-site partition semantics.
- Rotation/permutation invariance, numerical force/strain checks and custom radial
  compatibility. Combine center selection, ordinary pruning/masking and k range.
- Shared coefficients remain shared after fitting, saving/loading and AD. Do not
  compare a newly independently fitted decomposed model as if it were identical.

## Design matrices, priors and derivatives

- Selected/transformed A*c + baseline equals independent E/F/W under fit_pair true
  and false. Check Unitful labels and baseline subtraction after transformations.
- Synthetic small fits recover predictions; do not demand unique coefficients in
  rank-deficient cases. Installation refreshes inference caches and preserves masks.
- Check mapped prior action against R_old*T, including simple subset maps and a
  nontrivial mapped-coordinate example if supported. Unknown custom degree requires
  an explicit/documented prior policy, not an invented physical degree.
- Actual supported Zygote/ForwardDiff calls for E->X/ps, F->X/ps, W->ps and existing
  W->X behavior. Independent finite differences of scalar objectives with fixed
  force/virial cotangents; no need for full Hessians/Jacobians.
- Include nonlinear custom parameters and radial scales, not readout-only checks.
  Pullbacks survive workspace reuse and caller mutation under existing ownership
  contracts. Repeated pullbacks and zero/thunk cotangents remain correct.

## Custom radial implementations

- Values, first distance derivatives and needed second/mixed derivative products
  checked against an independent analytic example and numerical stencils.
- Scalar versus lane-wise SIMD values/derivatives, representative widths, Float32
  and Float64, zero channel scales, safe inactive lanes, cutoff interior/boundaries
  and declared origin policy. Never rely on zero times NaN masking.
- Analytic joint evaluator and scalar AD fallback agree. Clearly label fallbacks
  with allocations/performance limits; test allocating adapters separately.
- Incorrect output shape/order declaration, unsupported functions or derivative
  paths, parameter shapes and double-envelope mistakes fail or are diagnosed.
- Imported spline path retains reference predictions. The optional ACEpotentials
  adapter uses pinned actual methods and agrees on sampled channels/derivatives.
- CPU scalar/SIMD/AD capability results are separate. GPU requirements are documented;
  no CPU check may set a verified GPU-capable flag without actual device compilation.

## Persistence

- Old native/resolved formats still load; all new selection and radial routes
  round-trip after relocation and in a fresh process. Restore authoritative ps/st,
  not stale precontracted weights or workspace contents.
- Source callable with nested config and custom trainable ps reproduces E/F/W and
  checked derivative products. Namespace collisions do not alter unrelated modules.
- A source bundle with a visible sentinel side effect must NOT execute when code
  is disallowed. Opt-in loading, wrong entrypoint/interface/hash/reference values,
  malformed source and missing dependencies receive explicit diagnostics.
- Invalid saves leave an existing destination byte-identical. Unsupported callable
  saving does not silently use Julia closure serialization or lose parameters.
- Tabulation uses independent off-grid samples for value and derivative errors,
  correct domain/cutoff behavior and recorded tolerance/precision. Check relevant
  second-derivative limitations and explain that original nonlinear parameter
  gradients do not survive fixed-table export.

## Optional plotting, measurements and reporting

- Core load without UnicodePlots; plot construction after loading the extension;
  correct underlying counts/curve samples and empty/filtered behavior. Return plot
  objects. Avoid exact terminal glyph assertions and large dependencies.
- Measure prepared kernel allocations after warmup, model/scratch structural sizes
  and representative species-order overhead. Separate preparation/public ownership
  costs. Preserve one shared sorted force array; report retained scratch honestly.
- Run full supported-version suites where available, including threaded execution
  and optional extensions. Keep reference environment intact. Record versions,
  actual thread count, commands and any unavailable checks in STATUS.md.
- Completion report states implemented APIs, numerical evidence and material limits.
  No general speedup, arbitrary AD order or GPU support claims from these checks.
