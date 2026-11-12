#!/usr/bin/env bash
set -euo pipefail
# Run from repository root. Requires the supplied local data/ref_model inputs.
export JULIA_DEPOT_PATH=/tmp/pfw-may30-depot:$HOME/.julia
may30_julia=/var/home/teemu/.julia/juliaup/julia-1.13.3+0.x64.linux.gnu/bin/julia
may30_julia110=/var/home/teemu/.julia/juliaup/julia-1.10.12+0.x64.linux.gnu/bin/julia
may30_results=benchmark/cpu_ace/results/2026-11-11-may30
mkdir -p /tmp/pfw-may30-env
cp data/ace_ref/Project.toml data/ace_ref/Manifest.toml /tmp/pfw-may30-env/
JULIA_PKG_PRECOMPILE_AUTO=0 "$may30_julia" --startup-file=no --project=/tmp/pfw-may30-env benchmark/cpu_ace/may30_setup.jl
# Correctness and compatibility environments reuse the previous milestone's
# dependency snapshots, with the package path pointing at the current checkout.
for may30_version in 113 110; do
    mkdir -p "/tmp/pfw-may30-test${may30_version}-env"
    cp "benchmark/cpu_ace/results/2026-11-10-dag/environments/$may30_version/Project.toml" "benchmark/cpu_ace/results/2026-11-10-dag/environments/$may30_version/Manifest.toml" "/tmp/pfw-may30-test${may30_version}-env/"
done
"$may30_julia" --project=/tmp/pfw-may30-test113-env -t4 test/runtests.jl > "$may30_results/tests-113.log" 2>&1
"$may30_julia110" --project=/tmp/pfw-may30-test110-env -t4 test/runtests.jl > "$may30_results/tests-110.log" 2>&1
"$may30_julia" --project=/tmp/pfw-may30-test113-env -t4 benchmark/cpu_ace/may30_regression.jl > "$may30_results/focused-113.log" 2>&1
"$may30_julia110" --project=/tmp/pfw-may30-test110-env -t4 benchmark/cpu_ace/may30_regression.jl > "$may30_results/focused-110.log" 2>&1
# Optional base loading check: minimal environment, no angular/vendor package.
mkdir -p /tmp/pfw-may30-cpu-load-env
printf '[deps]\nPotentialFrameWork = "069ed237-161a-43fa-b820-4bdfd41b3db4"\n' > /tmp/pfw-may30-cpu-load-env/Project.toml
"$may30_julia" --project=/tmp/pfw-may30-cpu-load-env -e 'using Pkg; Pkg.develop(path=pwd())'
"$may30_julia" --project=/tmp/pfw-may30-cpu-load-env -t1 benchmark/cpu_ace/may30_cpu_load.jl > "$may30_results/cpu-only.log" 2>&1
# Only after correctness processes have finished: sequential warmed trials,
# one physical core, one Julia thread, BLAS forced to 1 by may30.jl.
taskset -c 0 "$may30_julia" --project=/tmp/pfw-may30-env -t1 benchmark/cpu_ace/may30.jl "$may30_results" > "$may30_results/benchmark.log" 2>&1
PFW_MAY30_RETIME=true PFW_MAY30_SECONDS=120 taskset -c 0 "$may30_julia" --project=/tmp/pfw-may30-env -t1 benchmark/cpu_ace/may30.jl "$may30_results" > "$may30_results/retime.log" 2>&1
python3 benchmark/cpu_ace/may30_report.py "$may30_results"
