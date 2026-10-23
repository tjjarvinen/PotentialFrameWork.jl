# Reference fixture and repository baseline

## Local reference inputs

Use `data/ace_ref/` in the user's local repository:
- `Project.toml`
- `Manifest.toml`
- `reference_model.json`
- `reference_results.xyz`
- `generate_reference_data.jl`

Read all of these. Do not assume that GitHub contains this directory. Preserve
the supplied files and model; generate additional exports in clearly named files.
Do not rerun the fitting script just to reload predictions: it fits and overwrites
the saved model. Load the existing model instead.

## Inspected fixture

The corrected archive was inspected statically on 22 October 2026:
- Julia 1.13.2; ACEpotentials 0.10.2; ACEfit 0.3.1.
- EquivariantTensors 0.4.3; SpheriCart 0.2.4; Polynomials4ML 0.5.10.
- Lux 1.31.4; LuxCore 1.5.3; AtomsBase 0.5.3;
  AtomsCalculators 0.2.3; AtomsSystems 0.1.2; ExtXYZ 0.2.2.
- JSON top-level keys: `manifest`, `meta`, `model_parameters`, `model_spec`,
  `project`. Embedded and external package versions agree.
- `model_spec`: model_name="ACE1", elements=["Ti","Al"], order=3,
  totaldegree=6, rcut=5.5 Å; Eref(Ti)=-1586.0195 eV,
  Eref(Al)=-105.5954 eV.
- `model_parameters` contains 270 numbers. This is a flattened parameter vector,
  not proof that all 270 entries are many-body invariant coefficients. Inspect
  reconstruction, component blocks and parameter ordering.
- XYZ: 33 structures, 32 with 2 atoms and one with 54 atoms; all periodic;
  23 have skew cells. Cutoff exceeds half the cell dimensions in the small cells:
  multiple images and nonzero self-images are essential, not optional edge cases.

These are file/source inspection results, not a numerical reproduction claim.
Recheck local inputs in case they have changed.

## JSON reconstruction

The file is a model recipe plus fitted parameters, not an explicit basis dump.
Pinned ACEpotentials source reconstructs via `make_model(D["model_spec"])` and
then `set_parameters!(model,D["model_parameters"])`.
Use the actual API, e.g. `ACEpotentials.load_model(path)` in this pinned version
returns the model and saved dictionary. Verify locally before coding adapters.

ACE1 compatibility construction also creates radial and pair bases and applies
version-specific defaults absent from the JSON. Resolve those defaults from the
constructed object. Export everything necessary to evaluate independently:
radial recurrences/transforms/envelopes, species channels, angular conventions,
product and coupling tables, numerical parameters, reference/pair components,
output scaling/transforms if any, and all parameter/basis index mappings.

Do not guess EquivariantTensors' actual involvement from its presence in Manifest.
Inspect which implementation the reconstructed object uses. Do not assume a
particular real/complex harmonic convention or radial l-dependence.

Prefer a versioned resolved-model export with explicit numerical tables. Keep
ACEpotentials construction in an optional adapter or dedicated reference/export
environment, with no repeated runtime dependency in our evaluator. Include source
model hash, producer versions and schema version in exported metadata. Reject
unsupported components with specific messages. No need to support every historic
ACEpotentials file format in this phase.

## Corrected force export and rounding

The original reference-generation script attached `forces` to system properties;
ExtXYZ retained the dataset's per-atom `force` values. The corrected script uses
`ExtXYZ.Atoms(sys)`, replaces `atom_data.force`, and updates system energy/virial.
Use the corrected version. The 33 frames' forces changed while the model,
environment, positions and cells remained identical. Stored net-force components
were about 1e-8 eV/Å or less, consistent with decimal rounding.

XYZ positions, cell and outputs have eight decimal places. Predictions were made
before geometry export; reloaded rounded geometries can therefore produce slightly
different results even in the original evaluator. For tight validation:
1. Load the saved model without refitting.
2. Read the saved XYZ geometries.
3. Recompute predictions on those exact geometries.
4. Save full-precision numerical references or compare live in the pinned env.
Do not impose machine-precision agreement with the rounded file or hide a genuine
mismatch by increasing tolerances indiscriminately. Keep dataset labels distinct
from model predictions and state units explicitly.

## Existing repository

Inspected remote commit: cd5d90e72b42834983151a9324ce84ed6d3b4029.
- NeighbourList owns storage/topology/scheduling, geometry, snapshot interface and
  neighbour benchmarks. PairPotentials owns models, evaluation, calculators and
  parameter energy. Root module includes/reexports the two functional modules.
- Full periodic-image and directed-snapshot support exists, including exact image
  labels relative to caller coordinates. Preserve the wrapping/winding logic.
- `geometry` enforces half-list ownership. Snapshot materialization adds both
  directions. Snapshot queries copy arrays; do not use them in final hot loops.
- `prepare!` rebuilds every call and allocates extraction/topology/schedule data.
  Some buffers are reused; there is no validated public Verlet/skin reuse.
- Existing pair force executor uses one sorted force array and phase scheduling.
- Explicit pair parameters currently support energy-only calculate requests.
  Generic parameter gradients use local ForwardDiff on scalar callbacks.
  SmoothPair's analytic gradient currently evaluates once per coefficient.
- Energy ChainRules rule explicitly excludes system derivatives and captures
  parameter gradients eagerly for safe workspace lifetime. It is not full AD.
- Production SmoothPair binds species and converts parameters to geometry precision;
  the separate scalar parameter path preserves promoted parameter arithmetic.
- AtomsCalculators outputs follow additive ! conventions and unitful boundary
  contracts. Preserve original-order outputs and existing custom pair callbacks.
- Latest STATUS entry records 14,970 assertions passing on Julia 1.10.12/1.13.2
  with four threads; this was not rerun for the instruction handoff.

Read current source/tests. Historical STATUS sections describe superseded geometry
limits and module boundaries; prefer current code and latest entries.
