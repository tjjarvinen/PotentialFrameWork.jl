# Exact stage-one reproduction

Run from the repository root. Application environments and their actual
Project/Manifest files are retained under `environments/`; adjust only their
local PotentialFrameWork source paths when reconstructing elsewhere. The
historical environment points at the immutable extracted baseline. Use the
original `harness-and-new-tests.tar.gz` for the primary recorded timing harness;
`final-harness-and-tests.tar.gz` additionally contains the controlled-repeat,
public-validation and report helpers and clarified legacy cache labels.

The original inputs remain untouched. Their hashes are in
`original-input-sha256.txt`. The baseline archive is reproducible from the
pinned Git commit as described in `../../baselines/README.md`.
The production code patch is `implementation.patch`, applied to that baseline.
Extract the harness/test archive over it to supply new tests and benchmark code.

```sh
export JULIA_DEPOT_PATH=/tmp/pfw-stage1-depot:$HOME/.julia
julia --project=/tmp/pfw-stage1-env -t4 test/runtests.jl
JULIA_DEPOT_PATH=/tmp/pfw-stage1-depot110:$HOME/.julia /var/home/teemu/.julia/juliaup/julia-1.10.12+0.x64.linux.gnu/bin/julia --project=/tmp/pfw-stage1-env110 -t4 test/runtests.jl
PFW_ROOT="$PWD" julia --project=/tmp/pfw-stage1-optional -t4 /tmp/pfw-stage1-optional-tests.jl

julia --project=/tmp/pfw-stage1-modern -t1 benchmark/cpu_ace/generate_modern.jl benchmark/cpu_ace/results/2026-11-10/models
PFW_CPU_SECONDS=60 julia --project=/tmp/pfw-stage1-modern -t1 benchmark/cpu_ace/oracle_modern.jl benchmark/cpu_ace/results/2026-11-10/models benchmark/cpu_ace/results/2026-11-10/oracles-final

python3 benchmark/gpu_ace/acquire.py --destination /tmp/pfw-stage1-legacy-inputs
# The pinned legacy environment needs Julia 1.10.12; Julia 1.13 lacks its bundled MbedTLS_jll.
export JULIA_DEPOT_PATH=/tmp/pfw-stage1-depot110:$HOME/.julia
/var/home/teemu/.julia/juliaup/julia-1.10.12+0.x64.linux.gnu/bin/julia --project=/tmp/pfw-stage1-legacy benchmark/gpu_ace/export_legacy.jl /tmp/pfw-stage1-legacy-inputs/Ti_combined_data_oct24_2023_ace_order_3_degree_20_rcut_6.0_cw_moreV.json benchmark/cpu_ace/results/2026-11-10/models/ti-resolved.toml
/var/home/teemu/.julia/juliaup/julia-1.10.12+0.x64.linux.gnu/bin/julia --project=/tmp/pfw-stage1-legacy benchmark/gpu_ace/export_legacy.jl /tmp/pfw-stage1-legacy-inputs/ti64_combined_data_Apr08_rcut_6.0_order_3_degree_15_custom_reg.json benchmark/cpu_ace/results/2026-11-10/models/tialv-resolved.toml
PFW_CPU_SECONDS=30 /var/home/teemu/.julia/juliaup/julia-1.10.12+0.x64.linux.gnu/bin/julia --project=/tmp/pfw-stage1-legacy -t1 benchmark/cpu_ace/legacy_oracle.jl /tmp/pfw-stage1-legacy-inputs/Ti_combined_data_oct24_2023_ace_order_3_degree_20_rcut_6.0_cw_moreV.json benchmark/cpu_ace/results/2026-11-10/oracles-final/ti-oracle.toml
PFW_CPU_SECONDS=30 /var/home/teemu/.julia/juliaup/julia-1.10.12+0.x64.linux.gnu/bin/julia --project=/tmp/pfw-stage1-legacy -t1 benchmark/cpu_ace/legacy_oracle.jl /tmp/pfw-stage1-legacy-inputs/ti64_combined_data_Apr08_rcut_6.0_order_3_degree_15_custom_reg.json benchmark/cpu_ace/results/2026-11-10/oracles-final/tialv-oracle.toml

# Restore the Julia 1.13 depot for CPU timings; run all expensive commands sequentially.
export JULIA_DEPOT_PATH=/tmp/pfw-stage1-depot:$HOME/.julia
julia --project=/tmp/pfw-stage1-env benchmark/cpu_ace/suite.jl --project /tmp/pfw-stage1-env --oracle-data benchmark/cpu_ace/results/2026-11-10/oracles-final --preset quick --samples 10 --seconds 30 --codegen true --profile true --baseline-project /tmp/pfw-stage1-baseline-env --output benchmark/cpu_ace/results/2026-11-10/final-quick
julia --project=/tmp/pfw-stage1-env benchmark/cpu_ace/suite.jl --project /tmp/pfw-stage1-env --oracle-data benchmark/cpu_ace/results/2026-11-10/oracles-final --preset full --sizes 8 --threads 1,4,16 --samples 10 --seconds 60 --features false --baseline-project /tmp/pfw-stage1-baseline-env --output benchmark/cpu_ace/results/2026-11-10/full-1024
julia --project=/tmp/pfw-stage1-env -t1 benchmark/cpu_ace/run.jl --preset full --models medium-ti,large-tial --sizes 18 --precisions Float64 --samples 10 --seconds 90 --diagnostics false --oracle-data benchmark/cpu_ace/results/2026-11-10/oracles-final --output benchmark/cpu_ace/results/2026-11-10/large-11664/current-t1.toml
julia --project=/tmp/pfw-stage1-env benchmark/cpu_ace/suite.jl --project /tmp/pfw-stage1-env --oracle-data benchmark/cpu_ace/results/2026-11-10/oracles-final --preset quick --models ti,tialv --sizes 4 --threads 1 --samples 10 --seconds 30 --features false --baseline-project /tmp/pfw-stage1-baseline-env --output benchmark/cpu_ace/results/2026-11-10/legacy-128
julia --project=/tmp/pfw-stage1-env -t1 benchmark/cpu_ace/run.jl --models small-tial --sizes 4 --widths 1,8 --samples 10 --seconds 30 --oracle-data benchmark/cpu_ace/results/2026-11-10/oracles-final --output benchmark/cpu_ace/results/2026-11-10/widths/current-t1.toml

# Controlled-affinity repeats on this host: CPUs 0–15 are distinct physical cores.
for t in 4 16; do
    if [ "$t" = 4 ]; then cpus=0-3; else cpus=0-15; fi
    taskset -c "$cpus" julia --project=/tmp/pfw-stage1-env -t"$t" benchmark/cpu_ace/pinned.jl current benchmark/cpu_ace/results/2026-11-10/models benchmark/cpu_ace/results/2026-11-10/oracles-final benchmark/cpu_ace/results/2026-11-10/pinned/current-t"$t".toml
    taskset -c "$cpus" julia --project=/tmp/pfw-stage1-baseline-env -t"$t" benchmark/cpu_ace/pinned.jl historical benchmark/cpu_ace/results/2026-11-10/models benchmark/cpu_ace/results/2026-11-10/oracles-final benchmark/cpu_ace/results/2026-11-10/pinned/historical-t"$t".toml
    PFW_CPU_SIZES=8 taskset -c "$cpus" julia --project=/tmp/pfw-stage1-modern -t"$t" benchmark/cpu_ace/oracle_modern.jl benchmark/cpu_ace/results/2026-11-10/models benchmark/cpu_ace/results/2026-11-10/oracles-pinned-t"$t"
done
julia --project=/tmp/pfw-stage1-env -t1 benchmark/cpu_ace/check_public.jl benchmark/cpu_ace/results/2026-11-10/models benchmark/cpu_ace/results/2026-11-10/oracles-final benchmark/cpu_ace/results/2026-11-10/validation/public.toml
python3 benchmark/cpu_ace/pinned_report.py benchmark/cpu_ace/results/2026-11-10
```

Choose a new output directory for repeat runs to preserve these recorded
observations. Replace cache paths consistently when doing so. Raw suite files
also retain exact Julia commands and complete process wall times; these wall
times include startup and are never used as prepared-evaluator timings.
The optional ACEfit test wrapper is archived with the validation logs. It runs
the existing selected-fit testset without requiring unrelated UnicodePlots.
