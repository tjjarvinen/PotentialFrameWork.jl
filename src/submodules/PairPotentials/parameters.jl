# NamedTuple trees preserve parameter names and scalar/array shapes for
# ForwardDiff and structured ChainRules tangents.
parameter_zero(x::Real)=zero(x)
parameter_zero(x::AbstractArray{<:Real})=zero(eltype(x))
parameter_zero(x::NamedTuple)=sum(parameter_zero,values(x);init=0)
parameter_zero(x)=throw(ArgumentError("parameter fields must be real scalars, arrays or nested NamedTuples"))
function parameter_seed(parameters::NamedTuple)
    isempty(parameters) && throw(ArgumentError("explicit parameters must be nonempty"))
    sum(parameter_zero,values(parameters))
end

parameter_flatten(x::Real)=[x]
parameter_flatten(x::AbstractArray{<:Real})=vec(copy(x))
parameter_flatten(x::NamedTuple)=vcat(map(parameter_flatten,values(x))...)
function parameter_rebuild(template,flat)
    offset=Ref(0)
    function rebuild(x)
        if x isa NamedTuple
            return NamedTuple{keys(x)}(map(rebuild,values(x)))
        elseif x isa Real
            offset[]+=1
            return flat[offset[]]
        else
            n=length(x); lo=offset[]+1; offset[]+=n
            return reshape(copy(flat[lo:offset[]]),size(x))
        end
    end
    rebuild(template)
end
parameter_scale(delta,ps::NamedTuple)=map(x -> parameter_scale(delta,x),ps)
parameter_scale(delta,x::Real)=delta*x
parameter_scale(delta,x::AbstractArray)=delta.*x

function AC.get_parameters(p::SmoothPair)
    n=length(p.labels)
    (;amplitudes=[p.amplitudes[i,j] for j=1:n for i=1:j])
end
function with_parameters(p::SmoothPair, parameters::NamedTuple)
    keys(parameters)==(:amplitudes,) || throw(ArgumentError("SmoothPair parameters must contain only amplitudes; cutoff is fixed"))
    a=parameters.amplitudes
    n=length(p.labels)
    a isa AbstractVector{<:Real} && length(a)==n*(n+1)÷2 ||
        throw(ArgumentError("amplitudes must be a real vector of length n*(n+1)/2"))
    all(isfinite,a) || throw(ArgumentError("amplitudes must be finite"))
    table=Matrix{eltype(a)}(undef,n,n)
    k=0
    for j=1:n, i=1:j
        k+=1
        table[i,j]=table[j,i]=a[k]
    end
    # Parameter precision and the fixed cutoff precision are independent.
    SmoothPair{eltype(a)}(copy(p.labels),table,p.radius)
end
AC.set_parameters!(p::SmoothPair, parameters::NamedTuple)=with_parameters(p,parameters)

@inline function cell_energy(w::Workspace{W}, bound, task, seed) where W
    energy=seed
    for i in w.ranges[task], image in w.images[task]
        r=w.ranges[image.cell]
        for j=max(first(r),i+(image.self ? 0 : 1)):W:last(r)
            _,block=geometry(w,i,j,last(r),image)
            for k=1:W
                lane(block.active,k) || continue
                # SIMD.jl vectors cannot store arbitrary dual numbers. Geometry
                # remains SIMD-packed; callbacks use scalar, parameter-promoted
                # values with zero geometry derivatives. No parameter conversion
                # to the workspace float type and no force storage are involved.
                b=PairBlock(block.central_species,(block.neighbour_species[k],),
                    lane(block.dx,k)+seed,lane(block.dy,k)+seed,lane(block.dz,k)+seed,true,true)
                energy+=bound(b,EnergyOnly()).energy
            end
        end
    end
    energy
end

"""`parameter_energy(sys, potential, parameters; kwargs...)` returns numeric eV.
This is the AD boundary for energy with respect to explicit numeric parameters.
System, species, physical cutoff, sorting and ownership are fixed with respect
to these parameters. Geometry/topology are nevertheless rebuilt on every call.
Callbacks run on scalar blocks for every workspace width; this is a correctness
path, not a promise of SIMD or AD performance. `parallel=true` uses cell phases.
"""
function parameter_energy(sys::AbstractSystem{3},p::AbstractPairPotential,parameters::NamedTuple;
                          neighbourlist=nothing,width=nothing,precision=nothing,parallel=false,
                          bruteforce_threshold=nothing,kwargs...)
    seed=parameter_seed(parameters)
    updated=with_parameters(p,parameters)
    cutoff(updated)==cutoff(p) || throw(ArgumentError("parameter updates must preserve the physical cutoff"))
    w=select_workspace(sys,neighbourlist,width,precision,bruteforce_threshold)
    prepare!(w,sys,cutoff(p))
    t=time_ns()
    bound=bind_species(updated,w.species_map)
    w.timings[:parameters]=(time_ns()-t)*1e-9
    t=time_ns()
    seed+=zero(eltype(w.x))
    partials=fill(seed,length(w.ranges))
    if parallel
        for phase in w.phases
            Threads.@threads :dynamic for k in eachindex(phase)
                task=phase[k]
                partials[task]=cell_energy(w,bound,task,seed)
            end
        end
    else
        for task in eachindex(w.ranges)
            partials[task]=cell_energy(w,bound,task,seed)
        end
    end
    w.timings[:evaluation]=(time_ns()-t)*1e-9
    t=time_ns()
    energy=sum(partials;init=seed)
    w.timings[:partial_merge]=(time_ns()-t)*1e-9
    w.timings[:output_permutation]=0.0
    energy
end

"""Gradient of numeric-eV `parameter_energy`, shaped like its parameter NamedTuple.
Extend for analytic parameter derivatives. Loading ForwardDiff supplies a generic
fallback; ChainRules uses this hook rather than tracing mutable neighbour storage.
"""
function energy_parameter_gradient(args...;kwargs...)
    throw(ArgumentError("load ForwardDiff for generic parameter gradients, or implement energy_parameter_gradient for this potential"))
end

# SmoothPair is linear in its independent symmetric amplitudes. Unit basis
# evaluations give exact derivatives without an AD dependency. This simple
# reference method rebuilds per coefficient; specialised fused gradients can
# replace it when large parameter tables become a target workload.
function energy_parameter_gradient(sys,p::SmoothPair,parameters::NamedTuple;kwargs...)
    with_parameters(p,parameters) # validate before allocating the basis
    a=parameters.amplitudes
    g=map(eachindex(a)) do k
        basis=zeros(eltype(a),length(a)); basis[k]=one(eltype(a))
        parameter_energy(sys,p,(;amplitudes=basis);kwargs...)
    end
    (;amplitudes=g)
end
