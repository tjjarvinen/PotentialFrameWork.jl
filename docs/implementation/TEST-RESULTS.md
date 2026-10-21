# Recorded validation — 2026-10-20

Hardware: AMD Ryzen 9 9950X, 16 cores / 32 logical CPUs, Linux x86_64.
Dependencies: AtomsBase 0.5.3, SIMD 3.8.1, Unitful 1.29.0. Julia 1.10 resolver
selected PrecompileTools 1.2.1 (the preinstalled 1.3 series was incompatible).
A writable temporary depot was layered over the existing read-only cache.

Final suite results (1,790 assertions in each process):

| Julia | Threads | Command | Result |
|---|---:|---|---|
| 1.10.12 | 1 | `julia --project test/runtests.jl` | pass |
| 1.10.12 | 4 | `julia --project -e 'using Pkg; Pkg.test()'` | pass |
| 1.10.12 | 32 | `julia --project test/runtests.jl` | pass |
| 1.13.2 | 4 | `julia --project test/runtests.jl` | pass |

Thread counts were set through JULIA_NUM_THREADS before starting separate Julia
processes; package-test subprocess inherited the setting. The 1.10.12 release
satisfies compat 1.10.10; the exact 1.10.10 patch release was not installed/tested.
Julia prerelease is left to the repository's existing CI matrix.

| Test set | Passing assertions |
|---|---:|
| Serial scalar oracle, species, units, permutation | 257 |
| Force and radial virial derivatives | 45 |
| Full Cartesian callback, nonsymmetric virial | 28 |
| SIMD widths, precisions, tails, cutoff, safe masking | 197 |
| Periodic ownership and actual write footprints | 655 |
| Parallel one-array execution and reset/reuse | 600 |
| Unsupported input errors | 8 |

Nonperiodic shear tests transform coordinates while retaining an inert diagonal
open cell: this is equivalent for open-boundary energy/forces, since the cell
has no image role. Public triclinic cell matrices are intentionally rejected.
Periodic dilation transforms both positions and cell. The nonsymmetric Cartesian
callback specifically tests the transpose/sign in the strain–virial relation.

`examples/custom_callback.jl` ran successfully on Julia 1.10.12 with four threads,
including both full-block and radial callbacks and energy-only evaluation.
`benchmark/inspect_simd.jl` ran successfully; LLVM contains `<4 x double>`
floating-point arithmetic. `git diff --check` passed. Benchmark output is stored
separately in `benchmark/results-1thread.txt` and `benchmark/results-4threads.txt`;
performance is deliberately not a CI assertion.


## ChemicalSpecies follow-up

The updated implementation passes 1,837 assertions (47 new species assertions)
on Julia 1.10.12 with 1 and 4 threads via Pkg.test, and Julia 1.13.2 with 4 threads.
These include explicit ChemicalSpecies field/callback types, isotope and atom-name
identity, species parameter lookup, permutation, and typed radial callbacks.
The custom example and SIMD LLVM inspection also passed. Earlier benchmark
records describe the pre-change integer-species implementation.

## AtomsCalculators integration

Final Pkg.test results: **2,332 passing assertions per process** on Julia 1.10.12
with 1 and 4 threads, and Julia 1.13.2 with 4 threads. AtomsCalculators 0.2.3 and
StaticArrays 1.9.23 are direct dependencies; other numerical dependency versions
are unchanged from the ChemicalSpecies follow-up.

New coverage totals 495 assertions: upstream AtomsCalculators conformance on
SmoothPair and RadialPotential (118), numerical/static-array/unit precision and
neighbour workspace creation/refresh checks (336), and additive mutating outputs,
fused tuple requests, positional/keyword state/parameters and energy-only calls
(41). Empty systems, explicit Float32 workspace overrides, changed positions,
species and cell, incompatible workspace options and destination sizes are tested.
The updated custom callback example (including AtomsCalculators calls) passed on
Julia 1.13.2 with four threads. No new performance claim is made.

## EmpiricalPotentials Lennard-Jones comparison

**2,616 assertions pass** in each full suite on Julia 1.10.12 with 1 and 4 threads,
and Julia 1.13.2 with 4 threads, against EmpiricalPotentials 0.2.6. The 284 new
checks compare energies, original-order forces and virials for Argon, using
AtomsSystems fixtures and the reference's inverse-distance cutoff tilt.
Scalar/SIMD widths, serial/threaded execution, cutoff edges, open/mixed/full
periodicity, small wrapped grids, force conservation, a finite difference,
immutability and singularity errors are covered. The 1.10 tests use
NeighbourLists 0.5.10 / AtomsCalculatorsUtilities 0.1.8; 1.13 uses 0.6.2 / 0.1.9.

The separate Julia 1.13.2 benchmark completed 224 cases spanning
256/864/4000/16384 atoms, 1/2/4/8 threads, scalar/SIMD and rebuild/fixed-topology
regimes. Every case passed energy/force/virial checks before and after timing.
Largest absolute discrepancies: 5.912e-12 eV energy, 4.025e-16 eV/Å force component,
5.059e-12 eV virial component. Raw records, reproduction instructions and plots
are in [the benchmark directory](../../benchmark/lennard_jones/README.md).

## Energy parameter AD follow-up, 2026-10-20

Full package tests: 2704 assertions passed on Julia 1.13.2 / 1 thread and
Julia 1.10.12 / 4 threads. The 88 new assertions cover explicit energy parameters,
ForwardDiff gradients/Hessians, analytic and finite-difference LJ parameter
references, ChainRules pullbacks and Zygote integration through calculate.
Existing species/periodic ownership, parallel footprints, force/virial and
EmpiricalPotentials comparison tests remain included. AD versions:
ForwardDiff 1.4.6, ChainRulesCore 1.26.4, Zygote 0.7.13.

The standalone example returned energy 0.0421875 eV and amplitude derivative
[0.421875] with both ForwardDiff and Zygote. The SmoothPair ChainRules rule also
passed in a fresh process without loading the ForwardDiff extension. The earlier
LJ timing matrix was not rerun and does not measure the new scalar AD path.
See STATUS.md for parameter contracts, implementation decisions and limitations.
