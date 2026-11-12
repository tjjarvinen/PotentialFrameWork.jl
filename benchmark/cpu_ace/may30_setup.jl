# Run with --project=/tmp/pfw-may30-env. Local Project/Manifest stay outside git.
using Pkg
Pkg.add(PackageSpec(name="ACEpotentials",version="0.10.3"))
Pkg.pin("ACEpotentials")
Pkg.develop(path=normpath(joinpath(@__DIR__,"../..")))
Pkg.add(["BenchmarkTools","SpheriCart","ACEfit","Lux","StaticArrays","AtomsBase","AtomsCalculators","ExtXYZ","Unitful"])
