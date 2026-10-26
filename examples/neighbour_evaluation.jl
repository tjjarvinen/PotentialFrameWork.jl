# julia --project examples/neighbour_evaluation.jl
using PotentialFrameWork, AtomsBase, Unitful, LinearAlgebra, StaticArrays
import AtomsCalculators as AC
const NI=PotentialFrameWork.NeighbourInterface
const NL=PotentialFrameWork.NeighbourList

# No mandatory potential supertype. The same callback accepts scalar or SIMD
# native pair blocks, and other providers implementing the accessor contract.
struct PolynomialPair
    radius::Float64
end
AC.length_unit(::PolynomialPair)=u"Å"
AC.energy_unit(::PolynomialPair)=u"eV"
NI.evaluation_kind(::PolynomialPair)=:pair
NI.search_cutoff(p::PolynomialPair)=p.radius
NI.species_channels(::PolynomialPair)=nothing
NI.prepare_potential(p::PolynomialPair,species;kwargs...)=p
NI.make_scratch(::PolynomialPair,capacity,precision)=nothing
function NI.evaluate_pair!(out,p::PolynomialPair,block,ps,scratch,props)
    r=NI.distance(block)
    q=1-r*r/p.radius^2
    any(p -> p isa AC.Energy,props) && NI.set_energy!(out,ps.a*q*q*q)
    if any(p -> p isa Union{AC.Forces,AC.Virial},props)
        NI.set_radial_derivative!(out,-6ps.a*q*q*r/p.radius^2)
    end
    nothing
end

# Squared vector density includes angular cross terms between distinct neighbors.
# Its small scratch is private to an execution slot; callbacks retain no views.
struct AngularDensity end
AC.length_unit(::AngularDensity)=u"Å"
AC.energy_unit(::AngularDensity)=u"eV"
NI.evaluation_kind(::AngularDensity)=:site
NI.search_cutoff(::AngularDensity)=3.0
NI.species_channels(::AngularDensity)=nothing
NI.prepare_potential(p::AngularDensity,species;kwargs...)=p
NI.make_scratch(::AngularDensity,capacity,::Type{T}) where T=zeros(T,3)
function NI.evaluate_site!(out,::AngularDensity,site,ps,scratch,props)
    fill!(scratch,0)
    for k=1:NI.neighbour_count(site)
        NI.active(site,k) || continue
        d=NI.displacement(site,k)
        for a=1:3; scratch[a]+=d[a]; end
    end
    any(p -> p isa AC.Energy,props) && NI.set_energy!(out,ps.a*sum(abs2,scratch))
    if any(p -> p isa Union{AC.Forces,AC.Virial},props)
        g=2ps.a*SVector{3}(scratch)
        for k=1:NI.neighbour_count(site)
            NI.active(site,k) && NI.set_gradient!(out,k,g)
        end
    end
    nothing
end
sys=AtomsBase.atomic_system([:H=>[0.,0.,0.]*u"Å",:He=>[1.,0.2,0.]*u"Å",:H=>[0.4,1.2,0.3]*u"Å"],
    ([10.,0.,0.]*u"Å",[0.,10.,0.]*u"Å",[0.,0.,10.]*u"Å"),(false,false,false))
ps=(;a=0.7)
pair=evaluate(sys,PolynomialPair(3.0);parameters=ps,request=(AC.Energy(),AC.Forces(),AC.Virial()))
input=NL.prepare_site_execution(AngularDensity(),sys;batch_size=2,width=4,slots=2)
NL.execute!(input.execution,input.positions,input.cell,ps;centers=[3,1])
angular=NL.site_result(input.execution)
@assert norm(sum(angular.forces;dims=2))<1e-12
println((;pair_energy=pair.energy,selected_angular_energy=angular.energy*u"eV"))
