# Compatibility entrypoints used by existing benchmarks and callers.
properties(::EnergyOnly)=(AC.Energy(),)
properties(::EnergyForcesVirial)=(AC.Energy(),AC.Forces(),AC.Virial())
properties(p)=NL.normalize_properties(p)
cell_task!(w,p,req,task)=NL.cell_task!(w,LegacyPair(p),properties(req),task)
execute!(w,p,req;kwargs...)=NL.execute!(w,LegacyPair(p),properties(req);kwargs...)

function outputs(w::Workspace{W,T},req) where {W,T}
    t=time_ns()
    energy=sum(w.energies)*u"eV"
    if req isa EnergyOnly
        w.timings[:partial_merge]=(time_ns()-t)*1e-9
        w.timings[:output_permutation]=0.0
        return (;energy)
    end
    v=reshape(vec(sum(w.virials;dims=2)),3,3)*u"eV"
    w.timings[:partial_merge]=(time_ns()-t)*1e-9
    t=time_ns()
    f=Matrix{typeof(zero(T)*u"eV/Å")}(undef,3,length(w.p))
    for k in eachindex(w.p), a=1:3
        f[a,w.p[k]]=w.force[a,k]*u"eV/Å"
    end
    w.timings[:output_permutation]=(time_ns()-t)*1e-9
    (;energy,forces=f,virial=v)
end

"""Evaluate an AtomsBase system, rebuilding all geometry/topology and returning
unitful original-order outputs. `evaluate!` reuses workspace storage, never the
caller's positions. Do not share one mutable workspace between concurrent calls.
"""
function evaluate!(w::Workspace,sys,p; request=EnergyForcesVirial(), parallel=false,parameters=nothing)
    props=properties(request)
    applicable(NI.evaluation_kind,p) && NI.evaluation_kind(p)!==:pair && throw(ArgumentError("pair evaluation requires a pair callback"))
    radius=applicable(NI.search_cutoff,p) ? NI.search_cutoff(p) : cutoff(p)
    prepare!(w,sys,radius)
    t=time_ns()
    legacy_parameters=!isnothing(parameters) && (p isa AbstractPairPotential || !applicable(NI.prepare_potential,p,w.species_map))
    legacy_parameters && (p=with_parameters(p,parameters))
    bound=if applicable(NI.prepare_potential,p,w.species_map)
        NI.prepare_potential(p,w.species_map;precision=eltype(w.x),width=NL.workspace_width(w))
    else
        LegacyPair(bind_species(p,w.species_map))
    end
    NL.check_units(bound)
    scratch=[NI.make_scratch(bound,NL.workspace_width(w),eltype(w.x)) for _ in w.ranges]
    w.timings[:parameters]=(time_ns()-t)*1e-9
    t=time_ns()
    NL.execute!(w,bound,props;parallel,parameters=legacy_parameters ? nothing : parameters,scratch)
    w.timings[:evaluation]=(time_ns()-t)*1e-9
    result=outputs(w,request isa Union{EnergyOnly,EnergyForcesVirial} ? request : props)
    result
end
function evaluate(sys,p; width=4, precision=Float64, bruteforce_threshold=32, kwargs...)
    evaluate!(Workspace(;width,precision,bruteforce_threshold),sys,p;kwargs...)
end

function outputs(w::Workspace{W,T},req::Tuple) where {W,T}
    props=NL.normalize_properties(req)
    t=time_ns()
    energy=NL.wants(AC.Energy,props) ? (;energy=sum(w.energies)*u"eV") : NamedTuple()
    virial=NL.wants(AC.Virial,props) ? (;virial=reshape(vec(sum(w.virials;dims=2)),3,3)*u"eV") : NamedTuple()
    w.timings[:partial_merge]=(time_ns()-t)*1e-9
    t=time_ns()
    forces=if NL.wants(AC.Forces,props)
        f=Matrix{typeof(zero(T)*u"eV/Å")}(undef,3,length(w.p))
        for k in eachindex(w.p),a=1:3; f[a,w.p[k]]=w.force[a,k]*u"eV/Å"; end
        (;forces=f)
    else
        NamedTuple()
    end
    w.timings[:output_permutation]=(time_ns()-t)*1e-9
    merge(energy,forces,virial)
end
