# Report formats (version 1)

All artifacts use TOML. Times are seconds unless a field explicitly says `ns`.
Numeric E/F/V errors use eV, eV/Å and eV respectively. Ratios must compare matching
model hashes, structures, precision and timed boundaries.

`PotentialFrameWork.GPUACESuite` stores options, subprocess commands, success
flags, log/report paths and omissions. A failed child causes a nonzero driver
exit after the other independent jobs finish. Paths identify the original run
location; distributed snapshots may relocate the associated files.

`PotentialFrameWork.GPUACEBenchmark` stores software/hardware versions, repository
commit/dirty status, source-file SHA-256 hashes, process/BLAS threads and elapsed
runtime. Each case contains model/export/structure hashes, species counts,
radial/angular/density/product/readout dimensions, precision, strategy, directed edge count, actual complete-call strategy, rounding
and absolute accuracy metrics. `prepared_EFV`, `complete_EFV` and optional
`trajectory` contain raw samples, median/minimum/maximum/standard deviation,
throughput, GC samples and **host** allocations. Device allocation counts are
explicitly unavailable. The nested `resident_array_payloads` contains the executor accounting categories.
Memory fields are array payloads, not a device allocator
high-water mark; shared payloads exclude compiler registers/spills and auxiliary
sorting/scan/skin allocations. `synchronized_stage_seconds` is an intrusive
single-call diagnostic, not additive evidence for undisturbed throughput.
`independent_oracle` records the original package's timing and exact boundary;
legacy SumIP uses three separate property calls. Missing entries are empty tables.

`PotentialFrameWork.CPUSIMDAudit` records hardware/compiler/package/source identity
and scalar/SIMD/angular/fusion/width cases with raw samples, host bytes, edge/lane
counts and colored phase sizes. Adjacent `.llvm.txt` and `.native.txt` artifacts
inspect warmed concrete methods. They are static evidence, not an estimate of the
fraction of execution time spent in SIMD. Compressed distributed artifacts retain
the original text when decompressed.

`PotentialFrameWork.PublicACEVerification` binds both resolved exports' byte
hashes to a successful independent CPU comparison. Standalone oracle fixtures
record source-model SHA-256, original package/Julia/thread versions, positions,
cell, species, E/F/V and raw timing samples. Public verification reports list
original-reference and rounded-CPU implementation errors separately. Derivative
reports identify the finite-difference step, device, tolerance and worst errors.
Raw test logs retain the assertions actually executed.

No schema assigns a universal Float32 suitability threshold. Implementation
checks and application precision sensitivity remain distinct.
