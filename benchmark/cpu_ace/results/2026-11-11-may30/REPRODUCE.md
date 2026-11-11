# May30 reproduction

Run [commands.sh](commands.sh) from the repository root. Use the Julia binaries
listed there or substitute equivalent installed versions. The local producer is
pinned to ACEpotentials **0.10.3**; it reloads the supplied checkpoint with
`load_may30_fit`, never fits or reads the original training data. All generated
outputs go to a distinct directory; the source checkpoint/script and historical
fixtures remain unchanged.

Producer and benchmark environment Project/Manifest files intentionally remain
under `/tmp/pfw-may30-env`, outside the deliverables. Package versions, source
file hashes and environment hashes are recorded with the results. The historical
core-test dependency snapshots referenced by the commands are separate from the
producer. Root Julia compatibility remains 1.10.10; a 1.10-compatible test snapshot
is used with Julia 1.10.12, rather than forcing modern ACEpotentials into it.

The original measured validation commands used `/tmp/pfw-may30-test-env` for the
1.13 snapshot and `/tmp/pfw-may30-test110-env` for 1.10. A fresh unchanged baseline
was obtained with `git archive fbbe8c4 | tar -x -C /tmp/pfw-may30-baseline`, copying
the 1.13 snapshot to `/tmp/pfw-may30-baseline-env` and replacing its PFW source path
with `/tmp/pfw-may30-baseline` in both Project and Manifest. Then:

```sh
JULIA_DEPOT_PATH=/tmp/pfw-may30-depot:$HOME/.julia /var/home/teemu/.julia/juliaup/julia-1.13.3+0.x64.linux.gnu/bin/julia --project=/tmp/pfw-may30-baseline-env -t4 /tmp/pfw-may30-baseline/test/runtests.jl
```

Live validation can be run without timings:

```sh
PFW_MAY30_TIMING=false PFW_MAY30_SIZES=4 JULIA_DEPOT_PATH=/tmp/pfw-may30-depot:$HOME/.julia /var/home/teemu/.julia/juliaup/julia-1.13.3+0.x64.linux.gnu/bin/julia --project=/tmp/pfw-may30-env -t1 benchmark/cpu_ace/may30.jl /tmp/pfw-may30-validation
```

The benchmark defaults to bcc replications `4,8,12` (128/1,024/3,456 atoms), plus
hcp and ω-like fixtures. `PFW_MAY30_SAMPLES=10`, `PFW_MAY30_SECONDS=10` bound each
trial. Actual sample counts and under-sampling flags must be inspected. The initial
ETACE 1,024/3,456-atom trials were under-sampled; a separate sequential run with
`PFW_MAY30_RETIME=true PFW_MAY30_SECONDS=120` repeats only those trials after
checking model hashes and polynomial parity again. Initial observations stay in
`initial-records.toml`; replacement trials stay in `records.toml`. Each median
uses one trial; samples from different trials are never pooled. The initial
harness snapshot is retained as `initial-harness.jl`. Alternative
sizes/node counts (`PFW_MAY30_SIZES`, `PFW_MAY30_NODES`) produce new representations
and must pass all parity checks before their performance is used. Use a new output
directory to preserve recorded measurements. Run expensive timed processes
sequentially, after correctness jobs finish. Only a serial comparison is claimed.

`data/may30/resolved.toml` is the complete exported model, and its three small
fixtures contain live splined ACEpotentials predictions. The result-directory
resolved-table link avoids storing a duplicate. Core regressions load through
`ACECalculator(path)` and do not load ACEpotentials/ETACE. The result fixtures
also retain the exact coordinates/cell/species used for every benchmark workload.
