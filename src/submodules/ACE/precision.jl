struct PrecisionReport{R,A,M}
    structures::R
    aggregate::A
    metadata::M
end
function Base.show(io::IO,::MIME"text/plain",r::PrecisionReport)
    println(io,"ACE precision comparison: ",length(r.structures)," structures; ",r.aggregate[:status])
    println(io,"  source ",r.metadata.source_precision,"; storage ",r.metadata.storage)
    for p in (:energy,:forces,:virial)
        haskey(r.aggregate,p) && println(io,"  ",p,": ",r.aggregate[p])
    end
end
function precision_evaluate(sys,b,props;kwargs...)
    if b.model isa GPUACEModel
        return gpu_result(execute_gpu!(GPUExecutor(sys,b.model,b.parameters;kwargs...);properties=props))
    end
    m=b.model; T=eltype(m.reference); inp=ace_inputs(precision_geometry(sys,T),m;precision=T)
    ex=CPUExecutor(m,inp.environment;precision=T,angular_backend=:block)
    execute_prepared!(ex,T.(inp.positions),T.(inp.cell),b.parameters;properties=props,kwargs...)
    result=prepared_result(ex)
    names=Tuple(p isa AC.Energy ? :energy : p isa AC.Forces ? :forces : :virial for p in props)
    NamedTuple{names}(Tuple(result[k] for k in names))
end

"""Compare Float32 and Float64 inference on one AtomsBase system or an iterable.
The same source coefficients/geometries produce both versions. Errors are absolute
in eV, eV/Å and eV. `tolerances=(energy=..., forces=..., virial=...)` optionally
checks total-energy difference and maximum component differences. No tolerances
means descriptive status only. Aggregate RMS pools components across structures;
failed/nonfinite cases make aggregate status failed and errors Inf. Empty inputs
return `:empty`. This is example sensitivity, not universal Float32 suitability.
"""
function compare_precision(m,ps,systems;array_type=Array,properties=(AC.Energy(),AC.Forces(),AC.Virial()),tolerances=NamedTuple(),kwargs...)
    props=NL.normalize_properties(properties)
    all(k in (:energy,:forces,:virial) && v isa Real && isfinite(v) && v>=0 for (k,v) in pairs(tolerances)) || throw(ArgumentError("tolerances must be finite nonnegative property errors"))
    systems=systems isa AtomsBase.AbstractSystem ? (systems,) : systems
    source=m isa GPUACEModel ? m.source : m
    sourceps=host_data(ps)
    bundles=map((Float64,Float32)) do T
        try
            (;bundle=adapt_model(source,sourceps;precision=T,array_type),error=nothing)
        catch e
            (;bundle=nothing,error=sprint(showerror,e))
        end
    end
    records=NamedTuple[]; totals=Dict{Symbol,Any}(); failed=false
    for (index,sys) in enumerate(systems)
        n=length(sys); results=Any[]; errors=String[]
        for (T,b) in zip((Float64,Float32),bundles)
            try
                isnothing(b.error) || error("preparation failed: "*b.error)
                result=precision_evaluate(sys,b.bundle,props;kwargs...)
                all(x->all(isfinite,x isa Number ? (x,) : x),values(result)) || error("nonfinite output")
                push!(results,result)
            catch e
                push!(results,nothing); push!(errors,string(T,": ",sprint(showerror,e)))
            end
        end
        bad=!isempty(errors); failed|=bad
        metrics=Dict{Symbol,Any}()
        for p in props
            key=p isa AC.Energy ? :energy : p isa AC.Forces ? :forces : :virial
            if bad
                metrics[key]=(;status=:failed,maximum=Inf,rms=Inf,component=nothing)
                totals[key]=(;maximum=Inf,squares=Inf,count=0,worst_structure=index,component=nothing)
                continue
            end
            x=results[1][key]; y=results[2][key]
            x=x isa Number ? [x] : vec(x); y=y isa Number ? [y] : vec(y)
            diff=Float64.(y).-Float64.(x); absdiff=abs.(diff)
            maxerr=maximum(absdiff;init=0.); component=isempty(diff) ? nothing : argmax(absdiff)
            rms=isempty(diff) ? 0. : sqrt(sum(abs2,diff)/length(diff))
            pass=!haskey(tolerances,key) || maxerr<=tolerances[key]
            status=haskey(tolerances,key) ? (pass ? :passed : :failed) : :descriptive
            failed|=!pass
            extra=key===:energy ? (;absolute=maxerr,per_atom=n==0 ? 0. : maxerr/n) :
                  key===:forces ? (;max_atom_vector=maximum(sqrt.(vec(sum(abs2,reshape(diff,3,:);dims=1)));init=0.)) : (;)
            location=key===:energy ? (;structure=index) : key===:forces ?
                (;structure=index,atom=isnothing(component) ? nothing : (component-1)÷3+1,axis=isnothing(component) ? nothing : (component-1)%3+1) :
                (;structure=index,row=isnothing(component) ? nothing : (component-1)%3+1,column=isnothing(component) ? nothing : (component-1)÷3+1)
            metrics[key]=merge((;status,maximum=maxerr,rms,component,location),extra)
            old=get(totals,key,(;maximum=-Inf,squares=0.,count=0,worst_structure=index,component=nothing))
            worse=maxerr>old.maximum
            totals[key]=(;maximum=max(old.maximum,maxerr),squares=old.squares+sum(abs2,diff),count=old.count+length(diff),
                worst_structure=worse ? index : old.worst_structure,component=worse ? component : old.component)
        end
        rounding=NaN; cellrounding=NaN
        identity=bytes2hex(SHA.sha256(string(typeof(sys),index)))
        try
            coords=[ustrip(u"Å",position(sys,i)[a]) for a=1:3,i=1:n]
            H=[ustrip(u"Å",cell_vectors(sys)[j][a]) for a=1:3,j=1:3]
            rounding=maximum(abs,Float64.(Float32.(coords)).-Float64.(coords);init=0.)
            cellrounding=maximum(abs,Float64.(Float32.(H[isfinite.(H)] )).-Float64.(H[isfinite.(H)]);init=0.)
            identity=bytes2hex(SHA.sha256(string(coords,H,[species(sys,i) for i=1:n])))
        catch e
            push!(errors,"geometry metadata: "*sprint(showerror,e)); bad=true; failed=true
            for key in keys(metrics)
                metrics[key]=(;status=:failed,maximum=Inf,rms=Inf,component=nothing)
                totals[key]=(;maximum=Inf,squares=Inf,count=0,worst_structure=index,component=nothing)
            end
        end
        push!(records,(;index,identity,atoms=n,status=bad || any(v.status===:failed for v in values(metrics)) ? :failed : isempty(tolerances) ? :descriptive : :passed,
            position_rounding_max_A=rounding,cell_rounding_max_A=cellrounding,errors,metrics))
    end
    aggregate=Dict{Symbol,Any}(:status=>isempty(records) ? :empty : failed ? :failed : isempty(tolerances) ? :descriptive : :passed)
    for (key,v) in totals
        aggregate[key]=(;maximum=v.maximum,pooled_rms=isfinite(v.squares) ? (v.count==0 ? 0. : sqrt(v.squares/v.count)) : Inf,
            components=v.count,worst_structure=v.worst_structure,component=v.component)
    end
    if haskey(aggregate,:energy)
        peratom=[get(r.metrics[:energy],:per_atom,Inf) for r in records]
        aggregate[:energy]=merge(aggregate[:energy],(;maximum_per_atom=maximum(peratom;init=0.)))
    end
    if haskey(aggregate,:forces)
        vectors=[get(r.metrics[:forces],:max_atom_vector,Inf) for r in records]
        aggregate[:forces]=merge(aggregate[:forces],(;maximum_atom_vector=maximum(vectors;init=0.),worst_vector_structure=isempty(vectors) ? nothing : argmax(vectors)))
    end
    PrecisionReport(records,aggregate,(;source_precision=string(source.coefficient_precision===Float32 || parameter_precision(sourceps)===Float32 ? Float32 : Float64),parameter_precision=string(parameter_precision(sourceps)),stored_precision=string(eltype(source.reference)),storage=string(array_type),
        model_sha256=source.source_sha256,properties=Tuple(string(typeof(p)) for p in props),tolerances,
        preparation_errors=Tuple(b.error for b in bundles),units=(;energy="eV",forces="eV/Å",virial="eV"),weighting="pooled components; failed cases propagate Inf"))
end
