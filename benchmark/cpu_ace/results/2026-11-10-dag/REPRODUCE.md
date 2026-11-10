# CPU product DAG reproduction

Starting commit: `4985499b6eb12421f512c3cb9d27671e5f76ec9b`.
The only initial untracked input was the supplied instruction document, copied
as `starting-instructions.md`. `starting-status.txt` additionally sees this newly
created result directory. Tracked source was clean. Immutable source archive:
`git archive 4985499b6eb12421f512c3cb9d27671e5f76ec9b -o /tmp/pfw-dag-stage1.tar`,
SHA256 `dc32268f71bd7a22454c146fd6f7fac74dfbe92265de97350510c4d6d9266d93`.
Extract in `/tmp/pfw-dag-stage1`; copy stage-one application environment and replace
only the PFW source path with that directory in Project and Manifest. Original
stage-one files, cached models/oracles and pinned manifests remain untouched.

Use direct Julia executables if juliaup cannot create its configuration lock in
the sandbox. No new dependency resolution is needed with the recorded application
environments. `commands.sh` is the exact sequential benchmark invocation script;
set `JULIA_DEPOT_PATH` as there. Source/model/geometry/environment identities,
raw observations and actual sample counts are retained with the results.

Fixed settings: width 4, SIMD radial on, fusion off, SpheriCart angular,
BLAS threads 1, perturbed bcc fixtures unchanged, CPU pinning via taskset. Ordinary
public combined calls reuse their calculator/plan but rebuild connectivity and
prepare the existing executor each call. `construct_and_public` additionally
constructs a fresh calculator, graph, scratch and owning outputs on every sample.
No global cache, no skin. `plan_build` isolates graph construction and
`executor_build` includes graph, angular basis, scratch and geometry preparation.

Benchmark defaults are ten warmed samples and a 120-second cap per trial;
`--samples`/`--seconds` override the budget. Every trial uses BenchmarkTools
`evals=1`, raw times, spread, host bytes and allocations; samples below ten are
explicitly flagged. Correctness/compilation warmup precedes production trials. `warmup_s` is the
measurement helper's first invocation and may already be compiled; it is not a
standalone compilation-time claim. The focused smoke results under /tmp are development checks only.

For full regression suites, use the same application environments:

```sh
JULIA_DEPOT_PATH=/tmp/pfw-stage1-depot:$HOME/.julia /var/home/teemu/.julia/juliaup/julia-1.13.3+0.x64.linux.gnu/bin/julia --project=/tmp/pfw-stage1-env -t4 test/runtests.jl
JULIA_DEPOT_PATH=/tmp/pfw-stage1-depot110:$HOME/.julia /var/home/teemu/.julia/juliaup/julia-1.10.12+0.x64.linux.gnu/bin/julia --project=/tmp/pfw-stage1-env110 -t4 test/runtests.jl
```

The structural compatibility check compares original product indices once per
prepared execution/update (and on executor re-preparation). Direct single-site
callbacks require immutable prepared model structure; topology changes require a
new prepared potential, as do density/channel changes. Coefficient changes refresh
weights without rebuilding the graph. Slots share plan arrays read-only. Legacy
array scratch constructors lazily grow numerical DAG buffers on first DAG use;
ordinary executor and generic preparation preallocate them before site loops.

`complete-validation.sh` runs the immutable-source full E/F/V export, current
three-algorithm comparison, concrete forward/reverse callee LLVM/native/warntype
inspection, and explicit arbitrary-order reassociation measurements. The export
covers 24 model/size/precision fixtures and the comparison 72 algorithm cases.
`dag_codegen_summary.py` summarizes static code evidence; `dag_report.py` derives
all observation tables and sample-count audits from raw TOML. Neither runs timings.

Generated code is preserved in `generated-code.tar.gz`: extract it in this result
directory to restore `codegen/` and `callee-codegen/`. The small
`callee-inspection.toml` and `codegen-summary.json` are directly readable without
extraction. Immutable outputs are `validation/stageone-outputs.toml.gz`; decompress
before rerunning only the comparison, or rerun the export script. Hashes in the
comparison refer to the decompressed original TOML bytes.

`source-and-harness.tar.gz` contains final src/ext/test trees, Project, current
benchmark/inspection/report helpers and documentation. `implementation.patch`
records tracked changes against the immutable stage-one commit; new files are
in the archive. Application manifests are in `environments/`; CPU-only load uses
the separately recorded root environment. Update local source paths only when
reconstructing. `artifact-sha256.txt` identifies retained artifacts. No old result,
model cache or reference environment is overwritten.
