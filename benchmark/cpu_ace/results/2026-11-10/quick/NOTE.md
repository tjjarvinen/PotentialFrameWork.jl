Exploratory quick run. Legacy dependency precompilation overlapped part of these
measurements. Final serial measurements are stored separately. The first code
inspection includes the product dispatcher (which calls both algorithms);
final code inspection additionally dumps the actual unrestricted/reference callees.
Initial diagnostic reduction includes clearing shared output arrays; final
microbenchmarks isolate the additive reduction itself. Source/harness hashes in
the raw reports identify the version that actually ran. Do not merge stage timings
across these revisions.
