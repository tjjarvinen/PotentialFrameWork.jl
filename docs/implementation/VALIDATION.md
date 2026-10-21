# Numerical, scheduling and performance validation

Tests should establish physical and ownership contracts, not merely reproduce
implementation steps. Use deterministic inputs/seeds where useful. Use an
independent all-pairs reference that does not share production cell traversal.

## Numerical behaviour

- Compare energy, forces and virial against the reference for several atom counts,
  species mixtures, uneven densities and cell occupancies.
- Compare forces with finite differences of energy away from singularities and
  cutoff discontinuities; use a smooth cutoff for this test.
- Check the chosen virial by a homogeneous strain finite difference: with
  x'=(I+epsilon)*x and H'=(I+epsilon)*H, the directional derivative is
  dE/dt = -sum(W .* transpose(epsilon)) under the specification's convention. Test at fixed
  neighbour topology away from cutoff crossings; test shear as well as dilation.
- Check total force is near zero for pair interactions and energy is translation
  invariant. Test permutation equivariance of forces and energy invariance.
- Test Float64 and Float32 with appropriate tolerances, widths 1 and representative
  hardware widths, and neighbour counts below/above block width.
- Test zero atoms if the system interface allows it, one atom, no neighbours,
  duplicate/self exclusions, exact cutoff policy and inactive lanes. Real coincident
  interacting atoms may produce a documented singularity; padding must not.
- Test energy-only and combined requests, repeated execution with reused outputs,
  and changed positions/species/cell according to the advertised invalidation rules.

## Boundary and ownership behaviour

- Nonperiodic, fully periodic and partially periodic orthorhombic cells.
- Atoms at faces and crossing periodic boundaries; small cell grids where stencil
  offsets wrap to the same cell; pairs near the minimum-image restriction.
- Every accepted unordered physical pair occurs once. No self-image interactions
  should be invented by wrapping. Test rejected cutoff/cell combinations explicitly.
- Brute-force check each phase's task write sets are disjoint on small cases.
  Assert actual destinations are subsets of advertised task footprints.
- Compare serial and threaded outputs repeatedly at several Julia thread counts,
  including more workers than tasks. Thread-count tests require separate Julia
  processes when changing startup thread count. Do not require bitwise equality
  unless deterministic accumulation order is explicitly implemented.

## Benchmarks

Provide scripts with reproducible commands. Warm up compilation, reset mutated
outputs between measurements, and prevent repeated force accumulation from
corrupting results. Keep timings separate from brittle CI performance assertions.

Cover small and large systems, uniform and inhomogeneous distributions, single
and multiple species, tail-heavy cell populations, multiple vector widths and
thread counts. State CPU, Julia/package versions, precision and affinity settings.

Measure complete public calls as well as stages: system extraction, permutation,
list construction, schedule construction, geometry packing/refresh, callback plus
force accumulation, energy/virial partial merge and final output permutation.
For reuse, state exactly which costs are amortised and over how many evaluations.
Record force storage separately from other metadata, peak/allocated memory where
available, masks/lane occupancy, phase count and task-work imbalance. Confirm that
force storage does not scale as atom count times thread count in the main path.

Use profiling/allocation inspection to guide optimisations. Inspect representative
SIMD kernel LLVM/native output for vector operations; do not demand architecture-
specific instruction names in portable tests. Record any spill or bandwidth concern
as measured evidence, not an assumed consequence of scalar versus vector code.
