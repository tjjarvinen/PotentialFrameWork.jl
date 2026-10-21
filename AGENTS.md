# PotentialFrameWork.jl — implementation guidance

This is Teemu Järvinen's working package name. Preserve the package name, UUID,
authorship, license and existing user changes. Develop components as submodules;
the neighbour component is `NeighbourList` (singular).

Read these files before implementation, in order:
1. `docs/implementation/INITIAL-SPEC.md`
2. `docs/implementation/DEVELOPMENT-PLAN.md`
3. `docs/implementation/VALIDATION.md`

Implement the initial finite-cutoff pair-potential stack described there. The
central goals are explicit SIMD pair blocks and spatially scheduled force
accumulation using one force array, with a small-buffer alternative later.
Sorting already works in NeighbourLists.jl; investigate and reuse its approach
rather than treating spatial sorting as an unproven research problem.

## Working rules

- Inspect the current repository, dependency APIs and existing instructions first.
  The handoff was based on a minimal package skeleton, not a completed design.
- Target the repository's Julia compatibility (currently `1.10.10`). Verify
  dependency compatibility and use actual package APIs, not invented methods.
- Use AtomsBase for system input. AtomsSystems is an optional helper, not a
  required concrete system type. Keep generic system access outside inner loops.
- Prefer SIMD.jl for the first CPU vector implementation. Keep width selectable
  through concrete types; do not promise GPU portability by widening CPU vectors.
- Keep atom/species data as runtime values, with concrete types in hot loops.
- Implement and validate serial execution before parallel scheduling.
- Make block evaluation and accumulation inlineable; avoid retaining all pair
  results in memory. Separate conceptual interfaces need not create runtime costs.
- Do not use a full force array per thread as the production implementation.
  Such an implementation is allowed as an explicitly labelled benchmark baseline.
- Do not implement ACE, MACE, training, GPU kernels, MPI or long-range electrostatics
  in this first version. Preserve extension points described in the specification.
- Do not spawn subagents unless the user explicitly requests delegation.
- Run meaningful numerical and scheduling tests. Report actual measurements,
  versions, hardware and thread counts; never invent speedups.
- Record material implementation choices and limitations in
  `docs/implementation/STATUS.md` as implementation proceeds. Do not silently
  weaken the specification to make a benchmark or test pass.

Proceed through the milestones without requesting approval for routine reversible
implementation decisions. If a scope change is necessary, explain it clearly.
