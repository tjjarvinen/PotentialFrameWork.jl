# PotentialFrameWork.jl — CPU ACE and differentiable pair implementation

Preserve package name, UUID, authorship, license, existing public behaviour and
unrelated user changes. This is Teemu Järvinen's working package name.

## Active scope and reading order

This file supersedes the first-phase prohibition on ACE and training. That
prohibition applied to the completed initial pair-potential implementation.
Existing first-phase numerical, ownership and one-force-array requirements remain.
Read, in order:
1. `docs/implementation/ACE-REFERENCE.md`
2. `docs/implementation/ACE-SPEC.md`
3. `docs/implementation/ACE-PLAN.md`
4. `docs/implementation/ACE-VALIDATION.md`
Then inspect existing source, tests, Project.toml and the latest entries of
`docs/implementation/STATUS.md`. Old INITIAL-SPEC/DEVELOPMENT-PLAN/VALIDATION files
are historical context; the ACE documents govern new scope where they differ.

Implement the CPU phase through M8 in ACE-PLAN.md. It includes importing the
provided ACEpotentials model, a separate ACE evaluator, SIMD radial evaluation,
preallocated numerical workspaces, complete species-grouped environments,
ACEfit integration, Lux support, and ChainRules derivatives for pair and ACE
models. Test the optional angular SIMD/fusion approaches after the baseline.

GPU implementation is a separate later project phase. Preserve backend/layout
boundaries without adding GPU dependencies, kernels or claiming GPU support.
MACE layers, MPI, long-range electrostatics and pressure-coupled differentiable MD
are also outside this phase.

## Working rules

- Reinspect the actual local tree: the user's `data/ace_ref` may be uncommitted.
  Never replace local files with the older remote checkout.
- Keep functional source components as submodules. Preserve `NeighbourList`
  (singular) and `PairPotentials`; add `ACE`. Keep the root module thin.
- Preserve Julia compatibility, currently 1.10.10, where practical. Verify new
  dependencies' actual APIs/compatibility. Keep the pinned Julia 1.13.2 reference
  environment separate. Do not raise core compatibility just to load a reference
  package or regenerate reference tables.
- Use AtomsBase as the system boundary; AtomsSystems is optional. Convert units
  explicitly at the boundary; numerical kernels use Å/eV conventions.
- Use explicit parameters/state and separate scratch ownership. Keep a convenient
  calculator wrapper and compatibility adapters for existing pair callbacks.
- Keep species/parameters as runtime data and hot loops concretely typed.
- Prefer SIMD.jl for CPU blocks. Validate a serial correctness path first.
- Half-list pair geometry and force schedules are not complete ACE environments.
  Build and test full directed neighbourhoods and their force-write footprints.
- Use one shared sorted force array in the production CPU executor. Full force
  arrays per worker are allowed only as labelled benchmark baselines. Small
  scratch buffers per exclusive execution slot are allowed.
- Standard sparse maps and in-place matrix multiplication are the baseline.
  TensorOperations.jl is optional, not a required core dependency. Bumper is
  optional for irregular CPU temporaries, not required for workspace ownership.
- Support AD into positions and explicit parameters; never silently detach them
  through stored displacements, precision conversions or stale inference caches.
- Required derivatives: energy -> positions/parameters; forces -> positions/
  parameters; virial -> parameters. Virial -> positions is optional and must
  explicitly report unsupported requests rather than return zero.
- Validate force pullbacks explicitly. An energy rrule alone is insufficient for
  force-loss training or differentiation through MD.
- Do not spawn subagents unless the user explicitly requests delegation.
- Run meaningful numerical, derivative, memory and scheduling checks. Record
  measurements and dependency/hardware versions; do not invent passing tests,
  performance gains or automatic-differentiation support.
- Update STATUS.md with completed gates, remaining limitations and reproduction
  commands. Keep historical measurements clearly labelled.

Proceed autonomously through routine reversible implementation choices and fix
failures within scope. If a dependency or missing input genuinely blocks progress,
finish independent work and report the precise blocker. Never weaken a numerical
contract silently to get a test or benchmark to pass.

## ACE interface phase — 24 October 2026

Implement the CPU interface phase in this reading order:
1. `docs/implementation/ACE-INTERFACE-SPEC.md`
2. `docs/implementation/ACE-INTERFACE-PLAN.md`
3. `docs/implementation/ACE-INTERFACE-VALIDATION.md`

These additions extend the completed CPU ACE phase, not replace its numerical,
unit, force-storage, scratch, AD lifetime or compatibility requirements. Preserve
existing supported model imports and native checkpoint loading.

Scope: descriptor inspection and optional UnicodePlots extension; masking and
structural pruning; restricted construction from an existing basis; central-index
selection; neighbor-species correlation-order projection; custom radial callables,
derivatives and persistence. Use original indices and explicit parameters.

No new site-energy API, changing regions, central-species selector, independent
multi-group species-order counters, general ACE basis generator, model-composition
framework, QM/MM physics, GPU implementation, or MACE layers in this phase.
Preserve boundaries that allow those future additions. Do not delegate to agents
without an explicit user request. Finish authorized reversible work autonomously.


## CPU neighborhood interface phase — 26 October 2026

This phase follows the completed ACE interface work. Read, in order:
1. docs/implementation/NEIGHBOUR-INTERFACE-SPEC.md
2. docs/implementation/NEIGHBOUR-INTERFACE-PLAN.md
3. docs/implementation/NEIGHBOUR-INTERFACE-VALIDATION.md

These instructions supersede earlier ownership/module-placement requirements only
where needed for a local pure protocol submodule and generic neighborhood
execution. Preserve earlier numerical, periodic-image, scratch, output ownership,
parameter authority, unit and AD contracts.

Implement only the CPU evaluation interface, pair/ACE adoption, single-site path,
batch fallback and associated validation. Keep existing colored reduction and
existing supported AD. Generic model state is not required by new callbacks;
retain state compatibility in existing calculator/Lux adapters.

AtomsCalculators remains unchanged and supplies its existing property types and
unit declarations. The experimental protocol contains declarations, documentation
and any necessary minimal immutable interface vocabulary, not execution or
allocation algorithms. Backend-specific preparation, packing, batch fallback,
scheduling and reduction stay in implementation modules.

No GPU implementation/dependencies, MACE layers, gather/atomic reduction,
optimized batched mathematics, optimized AD, skin reuse, different-cutoff
composition framework or AtomsCalculators proposal in this phase.
Do not delegate. Preserve existing public APIs through adapters, including legacy
pair callbacks and ACE constructors/imports/persistence.

### Neighborhood interface implementation record — 26 October 2026

The local experimental protocol is now `NeighbourInterface`; generic complete
storage/views, CPU loop batching, colored scheduling and force reduction live in
`NeighbourList`. Pair/ACE mathematics and compatibility adapters stay in their
model modules. `ACE.FullEnvironment` remains a compatibility binding. See
`docs/neighbour_evaluation.md` for final signatures, preparation and view lifetimes,
and `docs/implementation/STATUS.md` for validation/measurements. Preserve these
boundaries and the existing functional AD path in subsequent work. The numerical
prepared executor is not newly declared differentiable. This record does not
authorize any deferred GPU, reduction, AD-optimization or upstream-proposal work.

Pair protocol refinement: distance is the primary pair input, with radial dE/dr
outputs converted to Cartesian gradients in the backend. Displacement access and
Cartesian gradient outputs remain optional, preserving anisotropic and legacy
callbacks. Keep the functional pair AD adapters compatible.
