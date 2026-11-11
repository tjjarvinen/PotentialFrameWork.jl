# Fitted May30 ACEpotentials, ETACE and PFW CPU comparison

11 November 2026. The fitted Ti–Al–V checkpoint converts successfully to
polynomial ETACE and to PFW's resolved spline representation. All five generated
workloads pass the three independent numerical comparisons below. This is a
validated baseline; no subsequent performance optimization was implemented.

At 1,024 atoms, PFW prefix → DAG improves **134.144 → 109.380 ms public** and
**98.355 → 74.629 ms prepared**. ACEpotentials with the **same spline tables**
takes **102.047 ms public**, slightly faster than PFW DAG. Polynomial ETACE is
substantially slower on this fit. These results do not establish a universal
ranking across models, hardware or radial representations.

## Provenance and representation

Starting commit: `fbbe8c4efb96ad935d2c006b0b1971998b832748`, PR #16. Its
[product DAG report](DAG-REPORT.md) remains the historical performance baseline.
Tracked source started clean; only `data/ref_model/` and the supplied comparison
instructions were untracked. Both were preserved. No subagents were used.

The actual inputs are `data/ref_model/train_may30_ace.jl` and
`data/ref_model/may30_d10_n500.jls`. Checkpoint SHA256:
`31daf4c0ae8ad736fd8e156af9fe26116cabff6013e9ff2e4a107501188d5510`. Script SHA256:
`604f9e7499cd750b6229aa810f04918f285844e7bf3f71bdde77fa2e68d67e74`. `load_may30_fit` verifies ACEpotentials 0.10.3;
the checkpoint itself records Julia 1.13.3 / ACEpotentials 0.10.3. Its settings
are degree 10, order 3, cutoff 6 Å, 500 training frames, QR, smoothness 4,
and E/F/V weights 30/1/1. Ti/Al/V references are
−1593.6711788890334 / −111.43241493046436 / −2074.174048827439 eV.
No fitting or original training-data loading occurred.

The model has **2,155 invariants × 3 central species**, **156 density
components**, **5,652 products** (30 order-1, 746 order-2, 4,876 order-3),
114 radial and 10 pair channels, and 49 angular components through l=6.
The 2,155 × 5,652 coupling has 5,652 nonzeros; the DAG has 6,192 multiply nodes.
WB and Wpair shapes are 2,155 × 3 and 10 × 3: 6,495 explicit readout parameters.
Species channels remain **Ti, Al, V = 22, 13, 23** throughout conversion.
ETACE retains all three components: one-body/pair/many-body, with cutoffs 0/6/6 Å;
effective concrete types are recorded in `metadata.toml`.

The source uses **L2 real solid harmonics**, in `l*l+l+m+1` order. The modern
pair envelope is `PolyEnvelope1sR`, `((r/rc)^(-p)-1)*(1-r/rc)`, with p=1 here.
PFW now imports these conventions without changing the legacy spherical or ACE1
envelope defaults. Cartesian derivatives, scalar and SIMD radial paths,
SpheriCart/generated/block angular paths, functional evaluation, precision
conversion and native checkpoint persistence retain the convention. Unsupported
scalar-feature conversion and device packing fail explicitly.

Splines come exclusively from `Models.splinify(calc.model, calc.ps; nnodes=1000)`.
The fitted WB/Wpair and original reference object are retained; spline radial
parameter blocks are empty and spline-compatible Lux state is initialized.
The v0.10.3 calculator-level convenience route would reinitialize fitted
readouts, so it is deliberately not used. The exact source knots, transforms,
envelopes and piecewise cubic coefficients are exported. The **polynomial
source is unchanged**, and the final spline node count is **1,000**.

The complete, standalone [resolved model and small fixtures](../../data/may30/README.md)
are loaded by `ACECalculator(path)` without a live oracle dependency. The result
resolved-table link avoids duplicating the same 46 MB TOML. Its SHA256 is
`8adfb4b1583e1167d1186994a0420afb798fb8617d103143d640229ffe926dc0`. Producer Project/Manifest files stay outside committed
deliverables; their hashes, resolved versions and actual source-file hashes are
recorded with the results.

## Independent numerical comparisons

These are maxima over perturbed, strained/nonorthogonal periodic bcc
128/1,024/3,456-atom, hcp 150-atom and ω-like 240-atom cases, all in Float64.
Each has all three species, nonzero forces and nontrivial virials. An independent
27-image enumeration selects a V–V nearest neighbor of center 1, checked against
the prepared graph. The bcc-1,024 case exercises a nonzero image shift.
Face heights exceed twice the 6 Å cutoff; minimum distances exceed 2.60 Å.
Coordinates and perturbations are fixed analytic sin/cos expressions, with no
random global state. Exact cells, compositions, distances, neighbor counts,
geometry hashes and **per-fixture** errors are in the
[complete observations](results/2026-11-11-may30/OBSERVATIONS.md).

| Reference → candidate | Energy error eV/atom | Force component eV/Å | Virial component eV |
|---|---:|---:|---:|
| Polynomial ACE → polynomial ETACE | 1.364e-12 | 1.688e-14 | 8.731e-11 |
| Polynomial ACE → spline ACE | 9.095e-13 | 1.138e-08 | 4.532e-07 |
| Spline ACE → PFW (worst prefix/DAG/public/prepared) | 9.095e-13 | 3.194e-12 | 1.128e-10 |

Numerical-parity tolerances are 2.2e-10 eV/atom, 2e-9 + 2e-11·max|F|
eV/Å and 2e-8 + 2e-11·max|W| eV. Energy tolerance has no large-Eref relative
slack. Spline approximation has separate limits of 2.2e-6 eV/atom,
2e-5 + 2e-7·max|F| eV/Å and 2e-4 + 2e-7·max|W| eV. Nonfinite values and shape
mismatches fail. The measured approximation is far below these limits.

Angular values/Cartesian derivatives were compared to the live source basis at
multiple radii and directions. Source-spline radial values/derivatives were
checked for every species pair, including 5.9 Å, rc−1e-6, rc and rc+1e-6.
At exactly rc the source Interpolations Dual bounds check rejects a nonzero
transform tangent; both envelopes have analytically zero value and slope there,
which is checked separately. Near-cutoff derivatives use the live source.
`basis-validation.toml` retains the sample errors and this endpoint distinction
is not used to relax full E/F/V parity.

## Warmed public combined E/F/V

AMD Ryzen 9 9950X, CPU affinity **0** (one physical core), Julia **1.13.3**, one
Julia thread and BLAS=1. All results are CPU Float64. ACEpotentials 0.10.3,
EquivariantTensors 0.4.3, Polynomials4ML 0.5.10, SpheriCart 0.2.4, Lux 1.31.4,
ACEfit 0.3.1, AtomsBase 0.5.3, AtomsCalculators 0.2.3, SIMD 3.8.2 and
BenchmarkTools 1.8.0 were resolved. Full versions and CPU topology are recorded.

Median milliseconds, identical fitted parameters and geometry per row:

| Fixture | ACE polynomial | ETACE polynomial | ACE spline | PFW prefix spline | PFW DAG spline |
|---|---:|---:|---:|---:|---:|
| bcc-128 | 20.026 | 193.389 | 12.823 | 16.308 | 13.363 |
| bcc-1024 | 155.783 | 2166.301 | 102.047 | 134.144 | 109.380 |
| bcc-3456 | 520.796 | 7394.755 | 339.737 | 461.506 | 378.957 |
| hcp-150 | 18.011 | 144.312 | 10.876 | 14.934 | 11.406 |
| omega-240 | 38.809 | 350.454 | 21.738 | 30.221 | 24.494 |

All public boundaries are `AtomsCalculators.energy_forces_virial(sys, calc)`
and return owning unitful outputs. Polynomial/spline ACEpotentials build a fresh
PairList per call (default Folds executor with one task); ETACE builds a fresh
interaction graph per call, shared between pair and many-body components, with
no cross-call graph reuse. Its one-body component needs no neighbor graph.
PFW rebuilds full directed environments, species packing and colored schedules
per call, refreshes authoritative weights, reuses numerical capacity and returns
original-order owning outputs. Both PFW algorithms explicitly use width 4,
SIMD=true, SpheriCart, fusion=false, parallel=false, slots=1. Effective algorithms
are checked; `:auto` remains prefix/reverse.

PFW-versus-spline-ACE compares identical exported radial tables. ETACE polynomial
and PFW spline evaluate the same fit with **different radial representations**;
the approximation error is reported independently above. No spline-ETACE timing
or equivalence claim is made.

## Prepared prefix versus DAG and allocation boundaries

| Fixture | Prefix ms | DAG ms | Prefix/DAG |
|---|---:|---:|---:|
| bcc-128 | 12.341 | 9.314 | 1.33× |
| bcc-1024 | 98.355 | 74.629 | 1.32× |
| bcc-3456 | 328.064 | 248.250 | 1.32× |
| hcp-150 | 11.643 | 8.118 | 1.43× |
| omega-240 | 22.764 | 17.132 | 1.33× |

Prepared trials time `execute_prepared!` with an unchanged complete environment,
including weight refresh, numerical E/F/V and shared-force accumulation. They
exclude geometry extraction, topology/scheduling, final total-energy/virial sums,
force permutation/copies and unitful output construction in `prepared_result`
and the public boundary. Returned executor buffers are borrowed. Both algorithms
use identical prepared inputs and settings. These are diagnostic timings,
not a substitute for the public measurements.

Every prepared trial reports **32 host bytes / one allocation** at this measured
call boundary; no zero-allocation claim is made. At 1,024 atoms, host allocation
volume per public call is about 516.6 MB polynomial ACE, 1,652.3 MB polynomial
ETACE, 256.5 MB spline ACE and 46.5 MB PFW (prefix and DAG equal). These are
allocated byte volumes, not resident or peak memory. PFW retains its single
shared sorted force array; no per-worker full force arrays were added.

Reconstruction, conversion, splinify, import, calculator/executor preparation
and plan construction are outside warmed call trials. Initial elapsed times,
including first-use compilation where applicable, are 9.054 s reload,
1.132 s ET conversion, 1.619 s spline generation and 1.915 s PFW import.
The separately warmed DAG plan build is 0.413 ms / 2,233,976 allocated bytes.
Per-case setup times are retained; first-prefix construction includes compilation
and must not be compared as an algorithmic constructor speed ratio.

All **35 final trials reached ten samples** (350 selected observations).
Two original under-sampled trials are also retained: 37 executed trials and
357 measured observations in total. Each trial uses evals=1 after warmup, sequentially. The initial ten-second
budget yielded only 5 ETACE samples at 1,024 atoms and 2 at 3,456. Those raw
trials remain in `initial-records.toml`; only these two were independently
repeated with a 120-second limit, after rechecking hashes and polynomial parity.
The final table uses replacement trials, without pooling their samples with the
initial trials. Raw samples, GC times, min/max spread, allocations and actual
counts are retained. See `summary.json` for the final sample audit.

## Regression outcomes and scope

The unchanged starting tree freshly passed **56,658 assertions / 62 testsets**
on Julia 1.13.3, four threads. The post-change full suites passed
**57,528 assertions / 63 testsets each** on Julia **1.13.3 and 1.10.12**, four
threads. The final focused imported-model suite passed **872 assertions** on
both versions, including two additional explicit force-loss VJP checks added
after the full runs: position and fitted WB/Wpair directional derivatives.
It also covers energy-position AD, force finite differences, a nonsymmetric
virial strain difference with unchanged edge-image topology, scalar/SIMD values
and slopes, all three angular backends, both public product algorithms,
precision/convention retention and native save/reload parity.

Existing scalar/SIMD, legacy import, persistence, geometry, scheduling, parameter
authority and functional AD suites pass. A separate minimal environment passed
CPU-only loading and solid block E/F/V with SpheriCart/vendor extensions absent.
No new core/oracle dependency or Julia compatibility increase was introduced.
Original checkpoint/script, historical fixtures and root/reference manifests
pass their recorded SHA256 checks. All required CPU conversion gates passed;
multicore scaling, spline ETACE and device execution are unmeasured optional or
out-of-scope cases. The new representation rejects unsupported device packing.

[Reproduction commands and limitations](results/2026-11-11-may30/REPRODUCE.md)
link the exact source identities, raw numerical/timing records, original
under-sampled trials and test logs. This milestone stops at the validated report.
