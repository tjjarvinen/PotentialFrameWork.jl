"""Shifted Lennard-Jones radial mathematics; cutoff filtering belongs to backend.
Use `lennard_jones(;epsilon,sigma,radius)` to construct an ordinary RadialPotential.
The origin is singular and explicitly rejected by GPU result checks.
"""
struct LennardJonesRadial{T}
    epsilon::T
    sigma::T
    radius::T
end
@inline function (p::LennardJonesRadial)(r,central,neighbors,request)
    T=eltype(typeof(r)); e=T(p.epsilon); s=T(p.sigma); c=T(p.radius)
    ratio=s/r; ratio2=ratio*ratio; q=ratio2*ratio2*ratio2
    ratio_c=s/c; qc2=ratio_c*ratio_c; qc=qc2*qc2*qc2
    energy=T(4)*e*(q*q-q-qc*qc+qc)
    request isa EnergyOnly && return (;energy)
    (;energy,derivative=T(24)*e*(q-T(2)*q*q)/r)
end
function lennard_jones(;epsilon=0.0103,sigma=3.4,radius=8.5)
    e,s,c=promote(float(epsilon),float(sigma),float(radius))
    all(isfinite,(e,s,c)) && e>0 && s>0 && c>0 || throw(ArgumentError("LJ parameters must be positive and finite"))
    RadialPotential(LennardJonesRadial(e,s,c),c)
end
AC.get_parameters(p::RadialPotential{<:LennardJonesRadial})=(;epsilon=p.radial.epsilon,sigma=p.radial.sigma)
function with_parameters(p::RadialPotential{<:LennardJonesRadial},ps::NamedTuple)
    keys(ps)==(:epsilon,:sigma) || throw(ArgumentError("LJ parameters require epsilon and sigma; cutoff is fixed"))
    lennard_jones(;ps.epsilon,ps.sigma,radius=p.radius)
end
function prepare_gpu_radial(p::LennardJonesRadial,::Type{T},ps) where T
    e=isnothing(ps) ? p.epsilon : ps.epsilon
    s=isnothing(ps) ? p.sigma : ps.sigma
    all(isfinite,(e,s)) && e>0 && s>0 || throw(ArgumentError("LJ parameters must be positive and finite"))
    LennardJonesRadial(gpu_number(T,e),gpu_number(T,s),T(p.radius))
end
