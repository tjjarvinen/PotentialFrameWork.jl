# Portable data only: never serialize environments, executors or inference caches.
const CHECKPOINT_SCHEMA = "PotentialFrameWork.ACECalculator"
const CHECKPOINT_TYPES = Dict(string(T)=>T for T in
    (Float32,Float64,Int8,Int16,Int32,Int64,UInt8,UInt16,UInt32,Bool))

function checkpoint_encode(x::NamedTuple)
    Dict("kind"=>"namedtuple", "names"=>string.(collect(keys(x))),
         "items"=>[checkpoint_encode(v) for v in values(x)])
end
checkpoint_encode(x::Tuple)=Dict("kind"=>"tuple", "items"=>[checkpoint_encode(v) for v in x])
checkpoint_encode(::Nothing)=Dict("kind"=>"nothing")
checkpoint_encode(x::String)=Dict("kind"=>"string", "value"=>x)
function checkpoint_encode(x::Union{Number,AbstractArray})
    T=x isa AbstractArray ? eltype(x) : typeof(x)
    haskey(CHECKPOINT_TYPES,string(T)) || throw(ArgumentError("unsupported checkpoint numeric type $T; save ordinary numeric parameters, not AD tracer values"))
    all(isfinite,x isa AbstractArray ? x : (x,)) || throw(ArgumentError("nonfinite checkpoint values"))
    if x isa AbstractArray
        return Dict("kind"=>"array", "eltype"=>string(T), "shape"=>collect(size(x)), "values"=>collect(vec(x)))
    end
    Dict("kind"=>"scalar", "eltype"=>string(T), "value"=>x)
end
checkpoint_encode(x)=throw(ArgumentError("unsupported checkpoint value $(typeof(x)); use nested NamedTuples/tuples, numeric arrays/scalars, strings or nothing"))

function checkpoint_decode(d)
    kind=d["kind"]
    if kind=="namedtuple"
        names=Symbol.(d["names"]); items=d["items"]
        length(names)==length(items) && length(unique(names))==length(names) || throw(ArgumentError("invalid checkpoint field names"))
        return NamedTuple{Tuple(names)}(Tuple(checkpoint_decode(v) for v in items))
    elseif kind=="tuple"
        return Tuple(checkpoint_decode(v) for v in d["items"])
    elseif kind=="nothing"
        return nothing
    elseif kind=="string"
        return String(d["value"])
    elseif kind in ("array","scalar")
        T=get(CHECKPOINT_TYPES,d["eltype"],nothing)
        isnothing(T) && throw(ArgumentError("unsupported checkpoint numeric type"))
        if kind=="scalar"
            value=T(d["value"])
            isfinite(value) || throw(ArgumentError("nonfinite checkpoint value"))
            return value
        end
        shape=Tuple(Int.(d["shape"]))
        all(>=(0),shape) || throw(DimensionMismatch("negative checkpoint array dimension"))
        vals=T[v for v in d["values"]]
        prod(big.(shape))==length(vals) || throw(DimensionMismatch("checkpoint array shape"))
        all(isfinite,vals) || throw(ArgumentError("nonfinite checkpoint array"))
        return reshape(vals,shape)
    end
    throw(ArgumentError("unsupported checkpoint value kind $kind"))
end

function checkpoint_radial(channels::AbstractMatrix{<:RadialChannel},spec;parameters=NamedTuple())
    entries=map(CartesianIndices(channels)) do ij
        ch=channels[ij]; t=ch.transform
        tr=Dict(string(k)=>getproperty(t,k) for k in fieldnames(typeof(t)))
        tr["normalization_rin"]=t.rin
        keys=ch.envelope===:PolyEnvelope2sX ? ("x1","x2","p1","p2","s") :
             ch.envelope===:ACE1_PolyEnvelope1sR ? ("rcut","r0","p") :
             throw(ArgumentError("unsupported checkpoint radial envelope $(ch.envelope)"))
        Dict("centre"=>ij[1],"neighbour"=>ij[2],"transform"=>tr,
             "envelope_kind"=>string(ch.envelope),"envelope"=>Dict(zip(keys,ch.envelope_parameters)),
             "nodes"=>ch.nodes,"coefficient_shape"=>collect(size(ch.coefficients)),"coefficients"=>vec(ch.coefficients))
    end
    Dict("kind"=>"piecewise_cubic_agnesi","spec"=>[Dict("n"=>n,"l"=>l) for (n,l) in spec],"channels"=>vec(entries))
end
function checkpoint_model(m::ResolvedModel{Float64},ps=NamedTuple())
    rows,cols,vals=findnz(m.coupling)
    Dict("schema"=>"PotentialFrameWork.resolved_ACE","schema_version"=>1,
         "source_sha256"=>m.source_sha256,"species_atomic_numbers"=>m.atomic_numbers,
         "reference_eV"=>m.reference,
         "angular"=>Dict("kind"=>"real_spherical","normalization"=>"L2","ordering"=>"l*l+l+m+1","lmax"=>m.lmax),
         "radial"=>checkpoint_radial(m.radial,m.radial_spec;parameters=radial_parameters(ps,:custom_radial)),
         "pair_radial"=>checkpoint_radial(m.pair_radial,m.pair_spec;parameters=radial_parameters(ps,:custom_pair)),
         "interface_policy"=>checkpoint_policy(m.policy),
         "density_indices"=>collect.(m.density),"products"=>m.products,
         "coupling"=>Dict("shape"=>collect(size(m.coupling)),"rows"=>rows,"columns"=>cols,"values"=>vals))
end
checkpoint_model(m)=throw(ArgumentError("checkpoint requires a resolved Float64 ACE model"))

function checkpoint_read_model(d,ps;allow_code=false)
    ps isa NamedTuple && hasproperty(ps,:WB) && hasproperty(ps,:Wpair) ||
        throw(ArgumentError("checkpoint parameters require WB and Wpair"))
    nz=length(d["species_atomic_numbers"])
    m=d["coupling"]["shape"][1]; nq=length(d["pair_radial"]["spec"])
    ps.WB isa AbstractMatrix && size(ps.WB)==(m,nz) &&
        ps.Wpair isa AbstractMatrix && size(ps.Wpair)==(nq,nz) || throw(DimensionMismatch("checkpoint readout dimensions"))
    for (key,n) in ((:radial,length(d["radial"]["spec"])),(:pair,nq))
        if hasproperty(ps,key)
            block=getproperty(ps,key)
            hasproperty(block,:scale) && block.scale isa AbstractArray && size(block.scale)==(n,nz,nz) ||
                throw(DimensionMismatch("checkpoint $key scale dimensions"))
        end
    end
    # Reuse the resolved-schema validator; parameters have one authoritative copy
    # in the checkpoint tree, retaining their original element types there.
    modeldata=copy(d)
    modeldata["parameters"]=Dict("WB_shape"=>collect(size(ps.WB)),"Wpair_shape"=>collect(size(ps.Wpair)),
        "ordering"=>"vec(WB), vec(Wpair); Julia column-major","flat"=>vcat(vec(ps.WB),vec(ps.Wpair)))
    first(read_resolved(modeldata;allow_code,custom_parameters=ps))
end

"""
    ACECalculator(path::AbstractString; resolved_path=..., kwargs...)

Load a self-contained checkpoint written by `save_calculator`, a resolved ACE
TOML export, or a recipe JSON with its matching resolved export. `resolved_path`
is used only for recipe JSON. Execution keywords (`width`, `slots`, `precision`,
`angular_backend`) configure fresh scratch; execution caches are never loaded.
The default angular backend requires SpheriCart to be loaded.
Custom source checkpoints are refused unless `allow_code=true` explicitly trusts
the saved Julia definitions. Registered/table channels need no code opt-in.
"""
function ACECalculator(path::AbstractString;resolved_path=joinpath(dirname(path),"resolved_model_v1.toml"),allow_code=false,kwargs...)
    model,ps,st=if lowercase(splitext(path)[2])==".json" && startswith(lstrip(read(path,String)),"{")
        load_recipe(path;resolved_path)
    else
        d=TOML.parsefile(path)
        if get(d,"schema",nothing)==CHECKPOINT_SCHEMA
            get(d,"schema_version",nothing) in (1,2) || throw(ArgumentError("unsupported ACECalculator checkpoint version"))
            ps=checkpoint_decode(d["parameters"]); st=checkpoint_decode(d["state"])
            checkpoint_read_model(d["model"],ps;allow_code),ps,st
        else
            read_resolved(d)
        end
    end
    ACECalculator(model,ps,st;kwargs...)
end

"""
    save_calculator(path::AbstractString, calc::ACECalculator;
                    parameters=calc.parameters, state=calc.state)

Save a versioned, self-contained TOML checkpoint of the resolved basis, explicit
parameters and state. Pass `parameters=trained_ps` when training used explicit
AtomsCalculators parameters without installing them in `calc`. Saving does not
modify the calculator. Existing files are replaced only after a complete write.

Nested NamedTuples/tuples, ordinary numeric arrays/scalars, strings and `nothing`
are supported. Array shapes and Float32/Float64 types are retained. Unsupported
objects, AD tracers and nonfinite values are rejected. Geometry, scheduling,
workspaces and derived coefficients are rebuilt when loading the checkpoint.
This is PotentialFrameWork's format, not an ACEpotentials recipe export.
"""
function save_calculator(path::AbstractString,c::ACECalculator;parameters=c.parameters,state=c.state)
    # Encode/validate before opening the destination, including unsupported state.
    p=checkpoint_encode(parameters); s=checkpoint_encode(state); m=checkpoint_model(c.model,parameters)
    validate_checkpoint_model(c.model,parameters)
    if all(ch -> ch isa RadialChannel,c.model.radial) && all(ch -> ch isa RadialChannel,c.model.pair_radial)
        checkpoint_read_model(m,parameters) # preserve the original full validator
    end
    d=Dict("schema"=>CHECKPOINT_SCHEMA,"schema_version"=>2,"writer_julia"=>string(VERSION),
           "model"=>m,"parameters"=>p,"state"=>s)
    destination=abspath(path)
    tmp,io=mktemp(dirname(destination))
    try
        TOML.print(io,d;sorted=true)
        close(io)
        mv(tmp,destination;force=true)
    finally
        isopen(io) && close(io)
        isfile(tmp) && rm(tmp)
    end
    path
end
