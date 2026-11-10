#!/usr/bin/env bash
set -euo pipefail
# Run from repository root, after correctness suites finish. No concurrent jobs.
export JULIA_DEPOT_PATH=/tmp/pfw-stage1-depot:$HOME/.julia
dag_julia=/var/home/teemu/.julia/juliaup/julia-1.13.3+0.x64.linux.gnu/bin/julia
dag_results=benchmark/cpu_ace/results/2026-11-10-dag
# Main serial complete grid and isolated real-site diagnostics.
taskset -c 0 "$dag_julia" --project=/tmp/pfw-stage1-env -t1 benchmark/cpu_ace/dag.jl --sizes 4,8 --profile true --codegen true --output "$dag_results/serial.toml" > "$dag_results/serial.log" 2>&1
# Representative small and large cases at 11,664 atoms in both precisions.
taskset -c 0 "$dag_julia" --project=/tmp/pfw-stage1-env -t1 benchmark/cpu_ace/dag.jl --models small-tial,large-tial --sizes 18 --algorithms optimized,dag --diagnostics false --output "$dag_results/large.toml" > "$dag_results/large.log" 2>&1
# Identical mathematical/settings grid on immutable stage-one source.
taskset -c 0 "$dag_julia" --project=/tmp/pfw-dag-stage1-env -t1 benchmark/cpu_ace/dag.jl --models small-tial,large-tial --sizes 4,8 --precisions Float64 --algorithms optimized --stageone true --diagnostics false --output "$dag_results/immutable-stageone.toml" > "$dag_results/immutable-stageone.log" 2>&1
# Parallel scheduler unchanged; core 0..15 are distinct physical cores on host.
for dag_threads in 1 4 16; do
    case "$dag_threads" in
        1) dag_cpus=0 ;;
        4) dag_cpus=0-3 ;;
        16) dag_cpus=0-15 ;;
    esac
    taskset -c "$dag_cpus" "$dag_julia" --project=/tmp/pfw-stage1-env -t"$dag_threads" benchmark/cpu_ace/dag.jl --models small-tial,large-tial --sizes 8 --precisions Float64 --algorithms optimized,dag --parallel true --diagnostics false --output "$dag_results/scaling-t$dag_threads.toml" > "$dag_results/scaling-t$dag_threads.log" 2>&1
done
