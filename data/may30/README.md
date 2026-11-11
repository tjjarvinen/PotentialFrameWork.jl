# Fitted May30 Ti–Al–V import fixture

`resolved.toml` is the complete resolved spline export of the supplied
`data/ref_model/may30_d10_n500.jls`, reconstructed by `load_may30_fit` from
`train_may30_ace.jl` under ACEpotentials **0.10.3**, Julia **1.13.3**. No training
was performed, and the original training dataset is unnecessary for this import.
`provenance.toml` records checkpoint settings, package versions and source hashes.
The source checkpoint SHA256 is
`31daf4c0ae8ad736fd8e156af9fe26116cabff6013e9ff2e4a107501188d5510`.

This is the reduced degree-10, order-3, 500-frame fit with a 6 Å cutoff and
Ti/Al/V species order (atomic numbers 22/13/23). It contains 2,155 invariants,
156 density components, 5,652 products, 114 radial channels and 10 pair channels.
The angular basis is L2 real **solid** harmonics through l=6, ordered
`l*l+l+m+1`. The pair envelope is modern `PolyEnvelope1sR`, distinct from the
legacy ACE1 envelope. Readouts, references, species/coupling/index order and
original spline transforms/envelopes are retained. The 1,000 source spline
nodes are preserved; polynomial coefficients are exported in exact piecewise
cubic coordinates, with independently checked interior samples and derivatives.

The three fixture TOMLs contain exact deterministic coordinates, nonorthogonal
periodic cells, all three species and **splined ACEpotentials predictions** in
eV, eV/Å and eV. They are predictions, not original training labels. Generated
structures include asymmetric perturbations and an explicitly verified nearest
V–V image; the producer checks face heights > 12 Å for the common single-image
regime. These generated geometries are conversion/benchmark fixtures, not a claim
about physical validation of the original fit.

Load with `PotentialFrameWork.ACE.ACECalculator("data/may30/resolved.toml")`
after loading SpheriCart, or select `angular_backend=:block` without that optional
dependency. Both explicit prefix (`:optimized`) and DAG (`:dag`) algorithms work;
automatic selection remains unchanged. Native calculator checkpoints preserve
the angular convention. The scalar feature adapter and device packing reject
this solid/modern-envelope representation explicitly; no new GPU support is
provided here.

Producer and independent live checks:
`benchmark/cpu_ace/may30.jl`. Core import/finite-difference/derivative regressions:
`test/ace_may30.jl`. See the [comparison report](../../benchmark/cpu_ace/MAY30-REPORT.md)
and [provenance and raw results](../../benchmark/cpu_ace/results/2026-11-11-may30/REPRODUCE.md).
ACEpotentials and SpheriCart attribution remains in
`docs/implementation/licenses/`; neither oracle is a new core dependency.
