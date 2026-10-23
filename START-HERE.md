# CPU ACE implementation handoff — 22 October 2026

## Install

Copy this bundle's `AGENTS.md` to the repository root, merging unrelated local
instructions. Copy its `docs/implementation/ACE-*.md` files into the matching
directory. Leave the original implementation documents and STATUS.md in place;
these new documents explicitly define the next phase.

The reference files are already in the user's local `data/ace_ref` directory.
This bundle deliberately does not duplicate them. Use the corrected version with
model-predicted per-atom forces, corresponding to the uploaded `ace_ref(1).zip`.

No potential implementation is included in this bundle. No Julia tests were run
while preparing it. The source inspection baseline was remote main commit
`cd5d90e72b42834983151a9324ce84ed6d3b4029` (21 October 2026).
The local tree may have further changes and is authoritative.

## Suggested initial prompt for VS Code Codex

Read AGENTS.md and all ACE documents it references. Inspect the current source,
tests and local data/ace_ref before editing. Implement the CPU phase through M8
of ACE-PLAN.md in small reviewable steps, preserving the existing pair API and
one-shared-force-array execution requirement. Start by validating the saved
reference model in its pinned environment, establishing explicit position/
parameter interfaces and proving the required pair-potential derivative products.
Reconstruct and export the resolved ACE basis, implement full species-grouped
environments and an independent serial evaluator, then add ACEfit integration,
SIMD/preallocation/parallel execution and complete Lux/ChainRules support.
Benchmark angular block/fusion variants only after baseline correctness. Update
STATUS.md as gates pass. GPU implementation remains a separate later phase.
Do not stop for routine confirmation between milestones.

## Review checkpoints

1. Verified reference predictions; position/parameter contract and pair derivative
   checks; explicit model/state/workspace boundaries.
2. Imported model with intermediate, energy, force and virial agreement.
3. ACEfit assembly and fits; optimized CPU serial/parallel evaluator.
4. Required Lux/AD paths and short differentiable fixed-cell MD validation.
5. Angular SIMD/fusion experiment results and documented backend choices.

These checkpoints structure implementation/review; they are not mandatory
permission pauses. Do not implement GPU or MACE to complete this CPU handoff.
