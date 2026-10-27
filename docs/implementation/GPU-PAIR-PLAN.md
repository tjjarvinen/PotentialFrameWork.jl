# GPU phase 1 milestones

Proceed in order. Validate each milestone before optimizing it. Do not stop at a
plan or skeleton when implementation is authorized. Preserve unrelated changes.

## G0 — inspect and establish a working device

- Read current AGENTS.md, interface specification/documentation, relevant pair
  models, geometry, tests, Project.toml and latest STATUS.md. Record the base commit.
- Verify dependency APIs and Julia compatibility. Choose optional GPU loading;
  CPU imports must not initialize a GPU or require a driver.
- On the user's machine, record Julia/package versions, AMDGPU diagnostics, OS,
  driver/runtime and RX 6800 XT identity. Allocate device arrays and run a small
  KA kernel and AK sort/scan/reduction with explicit completion checks.
- If hardware is inaccessible, finish independent CPU-backed/kernel checks and
  scripts, then report the exact real-device gate as unverified. Never substitute
  a CPU backend run for a claim of GPU validation.

## G1 — prepared data and radial execution

- Add GPU workspace/configuration, reusable buffers, precision selection,
  species identity conversion/channel tables and current-parameter refresh.
- Add accessor-compatible device pair/output views. Compile a small radial
  callback through the actual NI contract. Validate sign, units, requests and origin.
- Temporarily consume CPU-prepared directed connectivity on GPU to isolate force
  evaluation. Clearly label this development/reference path.
- Deliver baseline central-force ownership, energy/virial reductions, original-order
  results and a device-numerical path. Validate both Float32/Float64 against CPU.

## G2 — complete fresh sorted-cell GPU path

- Build cell keys, records, sort, aligned arrays, occupied-cell ranges and lookup
  on GPU. Support sparse/open geometry without vacuum-sized allocation.
- Traverse current surrounding cell ranges directly, deduplicating periodic cell
  aliases and excluding self by identity. Preserve device data residency.
- Replace the G1 transferred-connectivity production path. A large structure must
  run without CPU neighbor construction or distance-matrix materialization.
- Exercise buffer reuse/growth and unsupported geometry errors. Retain CPU support.

## G3 — selectable candidate-thread tiling

- Implement bounded central/candidate tiles with central data reuse. Reduce each
  center's force correctly and document cross-workgroup ownership/partial storage.
- Exercise empty/tail tiles, dense cells, several candidate chunks and small systems.
- Compare identical physical results with baseline and CPU. Measure complete
  requested E/F/W computation, not only candidate distance throughput.
- Keep the experimental option even if slower; report measured limitations and
  tune only justified launch/tile choices. Do not promise a speedup.

## G4 — optional buffered connectivity

- Add explicit skin>0 mode and invalidation, fixed-order updates, reference
  displacement checking and conservative candidate storage. Skin=0 stays default.
- Implement supported input position/winding conventions explicitly. If validity
  cannot be guaranteed, rebuild conservatively instead of using stale candidates.
- Validate evolving geometries against rebuild-every-step CPU and GPU references.
- Measure list overhead versus rebuild savings. Do not tune skin automatically.

## G5 — validation, benchmarks and documentation

- Complete GPU-PAIR-VALIDATION.md and run existing CPU/interface/AD regressions.
- Generate deterministic structures using AtomsBuilder, convert to AtomsSystems,
  then benchmark. These are benchmark/test dependencies, not required core inputs.
- Provide copy-pastable commands and executable examples for both strategies,
  precision overrides, parameter updates, requested properties and optional skin.
- Record limitations and numerical tolerances. Update STATUS.md and test-result
  documentation with actual outcomes, including any outstanding hardware gates.
- Do not claim untested vendor support, passing tests or allocation/performance
  properties. Do not weaken tests to obtain a benchmark number.
