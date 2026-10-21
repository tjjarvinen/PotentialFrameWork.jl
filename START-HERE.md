# Codex implementation handoff

Copy AGENTS.md into the PotentialFrameWork.jl repository root and copy
docs/implementation/ into its docs/implementation/ directory. Merge with any
existing instructions rather than replacing them blindly. START-HERE.md is a
handoff note and need not be committed.

Repository inspected: https://github.com/tjjarvinen/PotentialFrameWork.jl
Baseline HEAD observed on 2026-10-20:
`1c58f1128a6298979d31a79389807e60537cce49`.
It contained a package skeleton, empty source module and placeholder tests;
Project.toml specified Julia 1.10.10 compatibility and package UUID
069ed237-161a-43fa-b820-4bdfd41b3db4. Reinspect the working tree before editing.

Suggested opening prompt in VS Code Codex:

> Read AGENTS.md and the implementation documents it references. Implement the
> initial finite-cutoff pair-potential stack through milestone 5. Start by inspecting
> the repository and dependencies, then proceed in the documented order. Preserve
> unrelated changes. Use the NeighbourList submodule, explicit SIMD blocks and
> cell-task scheduling into one shared spatially ordered force array. Validate
> numerical results, periodic ownership and parallel write footprints. Record
> decisions, completed checks, benchmarks and limitations in STATUS.md. ACE/MACE,
> GPU execution and long-range forces remain follow-up work.

These documents contain instructions only: no potential implementation has been
written or benchmarked as part of this handoff. Decisions such as starting with
orthorhombic cells and rebuilding lists by default are explicit first-version
defaults proposed here, not claims that every detail was previously agreed.
