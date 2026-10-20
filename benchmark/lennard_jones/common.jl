"""Shared benchmark model and deterministic AtomsSystems fixtures; no production dependency."""
module ArgonComparison
using PotentialFrameWork, AtomsBase, AtomsSystems, StaticArrays, Unitful
import AtomsCalculators as AC
import EmpiricalPotentials as EP
export ArgonLJ, argon_models, argon_fcc, argon_pair, numeric, errors, assert_agreement

"""Benchmark-only single-species LJ matching EmpiricalPotentials 0.2.6.
U=epsilon*[f(s)-f(sc)-(s-sc)*f'(sc)], f(s)=s^12-2s^6, s=rmin/r.
This is a tilt in INVERSE DISTANCE, not the conventional force shift in r.
Numeric parameters use eV/Å; active coincident atoms are singular.
"""
struct ArgonLJ{T} <: AbstractPairPotential
    epsilon::T
    rmin::T
    radius::T
    fcut::T
    dfcut::T
end
function ArgonLJ(epsilon,rmin,radius)
    all(x -> isfinite(x) && x>0,(epsilon,rmin,radius)) || throw(ArgumentError("positive finite LJ parameters required"))
    sc=rmin/radius
    ArgonLJ(epsilon,rmin,radius,sc^12-2sc^6,12*(sc^11-sc^5))
end
PotentialFrameWork.cutoff(p::ArgonLJ)=p.radius
function PotentialFrameWork.bind_species(p::ArgonLJ,labels)
    all(s -> isequal(s,ChemicalSpecies(:Ar)),labels) || throw(ArgumentError("benchmark LJ only supports Argon"))
    p
end
@inline function (p::ArgonLJ)(b::PairBlock{W,V},request) where {W,V}
    T=eltype(V)
    r2=b.dx*b.dx+b.dy*b.dy+b.dz*b.dz
    any(b.active & (r2==0)) && throw(DomainError(0,"coincident LJ atoms are singular"))
    invr=inv(sqrt(r2))
    s=T(p.rmin)*invr
    s2=s*s; s4=s2*s2; s6=s4*s2
    epsilon=T(p.epsilon)
    energy=epsilon*(s6*s6-2*s6-T(p.fcut)-(s-T(p.rmin/p.radius))*T(p.dfcut))
    request isa EnergyOnly && return (;energy)
    # ds/dr=-s/r and central force=(dU/dr)*d/r.
    duds=epsilon*(12*(s6*s4*s-s4*s)-T(p.dfcut))
    g=-duds*s*invr*invr
    (;energy,fx=g*b.dx,fy=g*b.dy,fz=g*b.dz)
end

"""Argon example parameters (not a fitted model claim): epsilon=0.010323 eV,
sigma=3.405 Å, rmin=2^(1/6)*sigma, cutoff=2.5*sigma. emins=-epsilon in EP.
"""
function argon_models(;epsilon=0.010323,sigma=3.405,radius=2.5*sigma)
    rmin=2^(1/6)*sigma
    ours=ArgonLJ(epsilon,rmin,radius)
    z=AtomsBase.atomic_number(ChemicalSpecies(:Ar))
    reference=EP.LennardJones(Dict((z,z)=>-epsilon*u"eV"),
                             Dict((z,z)=>rmin*u"Å"),radius*u"Å")
    ours,reference
end

"""Deterministically rattled FCC Ar: N=4*m^3, a=5.26 Å. No RNG/global state.
With periodic axes use m>=4 so all periodic lengths exceed twice the cutoff.
"""
function argon_fcc(m; periodic=(true,true,true),amplitude=0.04)
    m>0 || throw(ArgumentError("positive repeat count required"))
    a=5.26
    basis=((0.,0.,0.),(0.,0.5,0.5),(0.5,0.,0.5),(0.5,0.5,0.))
    positions=SVector{3,typeof(0.0u"Å")}[]
    for iz=0:m-1, iy=0:m-1, ix=0:m-1, b in basis
        i=length(positions)+1
        r=SVector(ntuple(k -> a*((ix,iy,iz)[k]+b[k])+amplitude*sin(i*(0.37+k)),3))
        wrapped=SVector(ntuple(k -> periodic[k] ? mod(r[k],m*a) : r[k],3))
        push!(positions,wrapped*u"Å")
    end
    L=m*a
    vectors=(SVector(L,0.,0.)*u"Å",SVector(0.,L,0.)*u"Å",SVector(0.,0.,L)*u"Å")
    AtomsSystems.generic_system(fill(ChemicalSpecies(:Ar),length(positions)),positions;
                               cell_vectors=vectors,periodicity=periodic)
end

function argon_pair(r;periodic=(false,false,false))
    AtomsSystems.generic_system([:Ar,:Ar],
        [SVector(0.,0.,0.)*u"Å",SVector(r,0.,0.)*u"Å"];
        cell_vectors=(SVector(30.,0.,0.)*u"Å",SVector(0.,30.,0.)*u"Å",SVector(0.,0.,30.)*u"Å"),
        periodicity=periodic)
end

function numeric(r)
    f=[ustrip(u"eV/Å",r.forces[i][a]) for a=1:3,i in eachindex(r.forces)]
    (;energy=ustrip(u"eV",r.energy),forces=f,virial=ustrip.(u"eV",r.virial))
end
function errors(a,b)
    x=numeric(a); y=numeric(b)
    (;energy=abs(x.energy-y.energy),force=maximum(abs,x.forces-y.forces;init=0.),
      virial=maximum(abs,x.virial-y.virial))
end
function assert_agreement(a,b;rtol=2e-10,atol=2e-10)
    x=numeric(a); y=numeric(b)
    for property in (:energy,:forces,:virial)
        isapprox(x[property],y[property];rtol,atol) || error("LJ $property disagreement: $(errors(a,b))")
    end
    errors(a,b)
end
end
