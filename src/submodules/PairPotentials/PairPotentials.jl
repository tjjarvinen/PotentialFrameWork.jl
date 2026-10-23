"""Pair models, force evaluation, calculator integration and parameter derivatives."""
module PairPotentials
using AtomsBase: ChemicalSpecies, AbstractSystem, position
using StaticArrays: SVector, SMatrix
using SIMD: Vec
using Unitful
import AtomsCalculators as AC
using ..NeighbourList: Workspace, PairBlock, pack, lane, prepare!, geometry, select
using ..NeighbourList: neighbour_list
import AtomsBase

export AbstractPairPotential, EnergyOnly, EnergyForcesVirial, SmoothPair,
    RadialPotential, cutoff, bind_species, with_parameters,
    evaluate, evaluate!, parameter_energy, energy_parameter_gradient

include("models.jl")
include("evaluation.jl")
include("calculators.jl")
include("parameters.jl")
include("numerical.jl")
end
