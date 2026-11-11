# Standalone focused core regression (use an environment with test extras).
using PotentialFrameWork, Test, Unitful, AtomsBase, SpheriCart, TOML, LinearAlgebra, StaticArrays
import AtomsCalculators as AC
include(joinpath(@__DIR__,"../../test/fixtures.jl"))
include(joinpath(@__DIR__,"../../test/ace_may30.jl"))
