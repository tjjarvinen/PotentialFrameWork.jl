using PotentialFrameWork, AtomsBase, Unitful
import AtomsCalculators as AC
using SIMD: Vec

# Full block callback: no traversal, masks, endpoint reduction or units in hot loops.
struct SoftSpring{T} <: AbstractPairPotential
    radius::T
    amplitude::T
end
PotentialFrameWork.cutoff(p::SoftSpring)=p.radius
@inline function (p::SoftSpring)(b::PairBlock, request)
    q=1-(b.dx*b.dx+b.dy*b.dy+b.dz*b.dz)/(p.radius*p.radius)
    energy=p.amplitude*q*q*q
    request isa EnergyOnly && return (;energy)
    g=-6*p.amplitude*q*q/(p.radius*p.radius)
    (;energy,fx=g*b.dx,fy=g*b.dy,fz=g*b.dz)
end

system=AtomsBase.atomic_system([:H=>[0.,0.,0.]*u"Å",:He=>[1.,0.2,0.]*u"Å"],
    ([8.,0.,0.]*u"Å",[0.,8.,0.]*u"Å",[0.,0.,8.]*u"Å"),(true,true,true))
p=SoftSpring(2.,0.7)
w=Workspace(width=4)
result=evaluate!(w,system,p;parallel=true)
println(result)
println(evaluate!(w,system,p;request=EnergyOnly()))

# Radial adapter: r is a scalar or explicit SIMD vector, depending on width.
struct SoftRadial{T}
    c::T
    a::T
end
@inline function (p::SoftRadial)(r,central_species,neighbour_species,request)
    q=1-r*r/(p.c*p.c)
    energy=p.a*q*q*q
    request isa EnergyOnly && return (;energy)
    (;energy,derivative=-6*p.a*r*q*q/(p.c*p.c))
end
radial=RadialPotential(SoftRadial(2.,0.7),2.)
@assert evaluate(system,radial).energy ≈ result.energy

# AtomsCalculators can allocate its neighbour storage or use the same workspace.
@assert AC.potential_energy(system,p) ≈ result.energy
calculator_result=AC.energy_forces_virial(system,p;neighbourlist=w,parallel=true)
@assert calculator_result.energy ≈ result.energy
println(calculator_result)
