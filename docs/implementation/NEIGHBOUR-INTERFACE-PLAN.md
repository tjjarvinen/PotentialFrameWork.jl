# Implementation plan

## M0 — Inspect and baseline

Read AGENTS.md, relevant previous ACE specification/validation, and the three new
documents. Inspect HEAD, dirty state, package versions, data/ace_ref, CPU/site/pair
paths, public calculators and AD extensions. Do not reset to the reviewed hash.
Run existing relevant regression tests and a small prepared pair/ACE timing and
allocation baseline when Julia/dependencies are available. Record exact commands
and hardware. Baseline failures are recorded and investigated, not hidden.

## M1 — Local pure protocol

Create NeighbourInterface with declarations/docstrings/accessor vocabulary.
Document required/optional capabilities and buffer ownership. Use existing AC
properties/units; keep request normalization, fallbacks and allocation outside
the protocol. Add a short third-party implementer example of a pair callback and
complete-site callback plus backend scratch preparation.

Gate: protocol has no model/backend imports, traversal, kernels or algorithms.

## M2 — Backend-neutral environments and preparation

Introduce lightweight native pair/site/batch views and output buffers implementing
the accessors. Separate generic environment topology/geometry/coloring from
ACE channel/cutoff preparation. Preserve existing FullEnvironment entry points.
Prepare concrete model-specific evaluator and exclusive scratch through hooks.
Record topology version/capacity/layout and reject stale prepared execution as
currently required. Preserve full periodic semantics.

Gate: a small analytic angular site callback consumes views without knowing
Workspace, FullEnvironment or CPUExecutor field names.

## M3 — Pair and ACE adoption

Refactor native pair execution through pair callbacks, preserving old callback,
request and parameter entry points through adapters.
Extract CPU ACE site mathematics and backend reduction. Keep all existing angular
options and model transformations. AtomsCalculators wrappers retain current
return shapes, additive ! behavior and optional state. Do not change persistence.

Gate: native pair/ACE public and prepared results match the independent references;
callbacks do not write global forces; generic backend has no ACE type checks.

## M4 — Batch fallback

Add batch_size with validation and default 1. Implement direct single-site
execution and larger batch fallback in implementation code. Reuse storage and
scratch, preserve colored footprints and selected-center semantics.
Test a separate callback implementing specialized batch dispatch only to verify
selection, not to optimize ACE contractions.

Gate: several batch sizes and SIMD widths agree serially/threaded, including tail,
empty selections, empty environments and uneven species counts.

## M5 — Conformance and compatibility

Add an independent tiny O(N²) reference backend/view implementation in test or
example code, restricted to a clearly declared geometry subset. It must construct
connectivity independently and implement accessors with different storage.
Run the same analytic pair and angular-site callbacks against it and native
execution. Exercise ACE if reasonably supported without broad adapter work.

This demonstrates independence without committing to a production adapter for
another package. A full NeighbourLists-provider adapter is deferred; retain and
test the existing optional query extension. If useful, adapt its returned arrays
only as a small test. Do not expand scope into a second optimized executor.

Gate: cross-implementation callbacks agree and unsupported geometry/layouts fail
explicitly. Third-party instructions require no edits to native source.

## M6 — Validate and document

Run NEIGHBOUR-INTERFACE-VALIDATION.md, compare prepared evaluation with baseline,
fix avoidable allocation/dispatch regressions, and document measured changes.
Explain any extra edge-buffer memory traffic introduced by callback separation;
do not assume zero overhead or claim improved speed.

Update README/examples and docs/implementation/STATUS.md with final signatures,
module ownership, completed gates, reproduction commands and limitations.
Update AGENTS.md to reflect completed phase without deleting historical rules.
Keep interface experimental and no upstream proposal.

Stop after this phase. GPU design is planned with the user after the actual
interface is complete; do not start kernels, MACE, AD optimization or reduction
experiments.

## Completion report

Report changed responsibilities, compatibility, tests actually run, measured
overhead/memory and unresolved limitations. Distinguish preserved AD coverage from
any unsupported direct differentiation of mutating executor APIs.
Do not present a future work proposal as completed implementation.

