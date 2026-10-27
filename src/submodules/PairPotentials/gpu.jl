# Accept plain inference numbers, rejecting dual/tracked numeric wrappers explicitly.
function gpu_number(::Type{T},x) where T
    x isa Union{AbstractFloat,Integer,Rational} || throw(ArgumentError("GPU differentiation/tracked numbers are unsupported"))
    y=T(x)
    isfinite(y) || throw(ArgumentError("GPU numeric input must be finite in selected precision"))
    y
end


"""Prepare a symmetric, stateless, device-compatible radial callback and parameters.
Third party models may extend this host hook, returning (prepared, parameters,
ordered ChemicalSpecies labels). Prepared callbacks implement NI.evaluate_pair!
using scalar distance and radial output; no CPU scratch or dictionaries are allowed.
"""
function prepare_gpu_potential(p,species,::Type{T},parameters) where T
    throw(ArgumentError("GPU radial preparation is unsupported for $(typeof(p)); legacy/asymmetric callbacks require explicit adaptation"))
end
struct GPUSmooth{T}
    radius::T
end
@inline function NI.evaluate_pair!(out,p::GPUSmooth,b,ps,scratch,props)
    r=NI.distance(b); T=typeof(r); q=one(T)-r*r/(p.radius*p.radius)
    a=@inbounds ps[b.ci,b.cj]
    NL.wants(AC.Energy,props) && NI.set_energy!(out,a*q*q*q)
    if NL.wants(AC.Forces,props) || NL.wants(AC.Virial,props)
        NI.set_radial_derivative!(out,-T(6)*a*q*q*r/(p.radius*p.radius))
    end
    nothing
end
function prepare_gpu_potential(p::SmoothPair,species,::Type{T},parameters) where T
    bound=isnothing(parameters) ? p : with_parameters(p,parameters)
    for s in species
        any(isequal(s),p.labels) || throw(ArgumentError("GPU model missing species $s"))
    end
    GPUSmooth(T(p.radius)),gpu_number.(T,bound.amplitudes),copy(p.labels)
end

gpu_inputs(sys,p,::Type{T},config) where T = NL.gpu_inputs(sys,NI.search_cutoff(p),T,config)

"""Evaluate GPU radial inference on an explicitly supplied KA backend.
`precision=Float64` is explicit and independent of input units. Default fresh cells
are rebuilt every call. Results own original-order host arrays with Å/eV units.
No GPU differentiation is provided. Use NL.gpu_execute! for device-owned outputs.
"""
function evaluate_gpu!(w::NL.GPUWorkspace{T},sys,p;parameters=nothing,
        request=EnergyForcesVirial()) where T
    input=prepare_gpu_execution(w,sys,p;parameters)
    y=execute_gpu!(input;request)
    NL.KA.synchronize(w.backend)
    result=map(Array,y)
    all(v->all(isfinite,v),values(result)) || throw(DomainError("singular origin or nonfinite GPU pair callback output"))
    merge(haskey(result,:energy) ? (;energy=only(result.energy)*u"eV") : (;),
        haskey(result,:forces) ? (;forces=result.forces*u"eV/Å") : (;),
        haskey(result,:virial) ? (;virial=result.virial*u"eV") : (;))
end
function evaluate_gpu(sys,p,backend;precision=Float64,config=NL.GPUConfig(),kwargs...)
    evaluate_gpu!(NL.GPUWorkspace(backend;precision,config),sys,p;kwargs...)
end

"""Explicit opt-in for symmetric device radial functions. Return a stateless isbits
callable in selected precision. Extend for custom radial callables; this certifies
symmetric directed semantics and requires a finite regular origin or rejection.
"""
function prepare_gpu_radial(radial,::Type{T},parameters) where T
    throw(ArgumentError("radial callable needs prepare_gpu_radial certification for GPU execution"))
end
struct GPURadial{R}
    radial::R
end
@inline function NI.evaluate_pair!(out,p::GPURadial,b,ps,scratch,props)
    derivatives=NL.wants(AC.Forces,props) || NL.wants(AC.Virial,props)
    v=p.radial(NI.distance(b),NI.central_species(b),NI.neighbour_species(b),
        derivatives ? EnergyForcesVirial() : EnergyOnly())
    NL.wants(AC.Energy,props) && NI.set_energy!(out,v.energy)
    derivatives && NI.set_radial_derivative!(out,v.derivative)
    nothing
end
function prepare_gpu_potential(p::RadialPotential,species,::Type{T},parameters) where T
    radial=prepare_gpu_radial(p.radial,T,parameters)
    isbitstype(typeof(radial)) || throw(ArgumentError("GPU radial callable must be isbits"))
    GPURadial(radial),nothing,collect(species)
end


"""Prepare device inputs and reusable workspace; no host neighbor topology.
The returned input borrows reusable upload buffers and may be updated in-place through `input.positions` in
original atom order. Geometry must remain finite and cell/species unchanged;
reprepare after cell/species/cutoff changes. Numerical execution returns owning
device results, without force readback. No GPU AD is implemented.
"""
function prepare_gpu_execution(w::NL.GPUWorkspace{T},sys,p;parameters=nothing) where T
    NL.check_units(p)
    input=gpu_inputs(sys,p,T,w.config)
    prepared,ps,labels=prepare_gpu_potential(p,input.species,T,parameters)
    metadata=(input.lengths,input.periodic,input.counts,input.radius,input.identities,labels)
    if !isequal(w.metadata,metadata)
        w.cache=nothing
        w.metadata=metadata
    end
    w.lengths=input.lengths; w.periodic=input.periodic; w.counts=input.counts; w.radius=input.radius
    host=(input.positions,input.identities,collect(reinterpret(UInt64,labels)),ps)
    signature=map(a->isnothing(a) ? nothing : (eltype(a),size(a)),host)
    if isnothing(w.transfer) || w.transfer.signature!=signature
        arrays=map(host) do a
            isnothing(a) ? nothing : NL.KA.allocate(w.backend,eltype(a),size(a)...)
        end
        w.transfer=(;signature,arrays)
    end
    for (dst,src) in zip(w.transfer.arrays,host)
        isnothing(src) || copyto!(dst,src)
    end
    positions,identities,device_labels,parameters_device=w.transfer.arrays
    w.generation+=UInt(1)
    (;workspace=w,generation=w.generation,positions,identities,labels=device_labels,prepared,prepared_state=Ref(prepared),
      parameters=parameters_device,model=p,species=input.species)
end
function prepare_gpu_execution(sys,p,backend;precision=Float64,config=NL.GPUConfig(),kwargs...)
    prepare_gpu_execution(NL.GPUWorkspace(backend;precision,config),sys,p;kwargs...)
end
function execute_gpu!(input;request=EnergyForcesVirial(),parameters=nothing)
    w=input.workspace
    input.generation==w.generation || throw(ArgumentError("GPU input invalidated by a later preparation on this workspace"))
    prepared=input.prepared_state[]
    if !isnothing(parameters)
        prepared,ps,labels=prepare_gpu_potential(input.model,input.species,eltype(w.x),parameters)
        ps===nothing || copyto!(input.parameters,ps)
        input.prepared_state[]=prepared
    end
    NL.gpu_execute!(w,input.positions,input.identities,input.labels,prepared,input.parameters;
        properties=properties(request))
end
