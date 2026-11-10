#!/usr/bin/env bash
set -euo pipefail
# Run after commands.sh, from the repository root.
export JULIA_DEPOT_PATH=/tmp/pfw-stage1-depot:$HOME/.julia
dag_julia=/var/home/teemu/.julia/juliaup/julia-1.13.3+0.x64.linux.gnu/bin/julia
dag_results=benchmark/cpu_ace/results/2026-11-10-dag
taskset -c 0 "$dag_julia" --project=/tmp/pfw-dag-stage1-env -t1 benchmark/cpu_ace/dag_validate_stageone.jl export "$dag_results/validation/stageone-outputs.toml" > "$dag_results/validation/stageone-export.log" 2>&1
taskset -c 0 "$dag_julia" --project=/tmp/pfw-stage1-env -t1 benchmark/cpu_ace/dag_validate_stageone.jl compare "$dag_results/validation/stageone-outputs.toml" > "$dag_results/validation/stageone-comparison.log" 2>&1
taskset -c 0 "$dag_julia" --project=/tmp/pfw-stage1-env -t1 benchmark/cpu_ace/dag_codegen.jl "$dag_results/callee-codegen" > "$dag_results/callee-codegen.log" 2>&1
python3 benchmark/cpu_ace/dag_report.py "$dag_results"
taskset -c 0 "$dag_julia" --project=/tmp/pfw-stage1-env -t1 benchmark/cpu_ace/dag_reassociation.jl "$dag_results/reassociation.toml" > "$dag_results/reassociation.log" 2>&1
python3 benchmark/cpu_ace/dag_codegen_summary.py "$dag_results"
