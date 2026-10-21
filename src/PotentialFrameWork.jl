module PotentialFrameWork
include("submodules/NeighbourList/NeighbourList.jl")
include("submodules/PairPotentials/PairPotentials.jl")

using .NeighbourList: Workspace, PairBlock, validate_schedule, neighbour_list, neighbours,
    benchmark_neighbourlist
using .PairPotentials: AbstractPairPotential, SmoothPair, RadialPotential, EnergyOnly, EnergyForcesVirial,
    cutoff, bind_species, with_parameters, parameter_energy, energy_parameter_gradient, evaluate, evaluate!

export NeighbourList, PairPotentials, Workspace, PairBlock,
    AbstractPairPotential, SmoothPair, RadialPotential, EnergyOnly, EnergyForcesVirial,
    cutoff, bind_species, with_parameters, parameter_energy, energy_parameter_gradient, evaluate, evaluate!, validate_schedule,
    benchmark_neighbourlist, neighbour_list, neighbours
end
