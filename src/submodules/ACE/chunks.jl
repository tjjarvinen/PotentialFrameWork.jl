"""Exclusive mutable CPU inference state for serial site-chunk evaluation.

Construct with `ACEEvaluationState(calc)`, and pass as the fifth argument to
`AtomsCalculators.calculate`. Distinct states can share model/read-only parameters.
The model state is retained in `model_state`. The numerical mutation path adds
no AD rules; functional `ace_evaluate` and its supported pullbacks are unchanged.
"""
mutable struct ACEEvaluationState{E,P,S}
    environment::E
    prepared::P
    model_state::S
    chunks::Any
    scratch::Any
    scatter_capacity::Int
    batch_size::Int
    active::Bool
end
function ACEEvaluationState(c::ACECalculator;state=deepcopy(c.state),scatter_capacity=4096)
    scatter_capacity isa Integer && scatter_capacity>0 || throw(ArgumentError("scatter_capacity must be positive"))
    # The angular provider, derived weights and product workspace are state-owned.
    original=c.executor.prepared
    p=typeof(original)(c.model,deepcopy(original.basis),copy(original.weights),original.simd,original.fused,original.product_algorithm,original.product_plan)
    ACEEvaluationState(FullEnvironment(;width=c.environment.width),p,state,nothing,nothing,Int(scatter_capacity),c.executor.batch_size,false)
end
mutable struct ChunkACEScratch{P,PS,S,O,R}
    prepared::P
    parameters::PS
    numerical::S
    outputs::O
    properties::R
    ordinals::Vector{Int}
    centers::Vector{Int}
end
"""Build exclusive ACE scratch for a prepared evaluation state. Parameters remain
explicit and derived weights are refreshed before every pass."""
function create_scratch(c::ACECalculator,state::ACEEvaluationState;parameters=c.parameters,properties=(AC.Energy(),AC.Forces(),AC.Virial()))
    T=eltype(state.prepared.weights)
    capacity=maximum((state.environment.offsets[end,i]-state.environment.offsets[1,i] for i in eachindex(state.environment.species));init=0)
    if isnothing(state.scratch) || length(state.scratch.numerical.Rs)<capacity
        numerical=NI.make_scratch(state.prepared,capacity,T)
        outputs=[NL.SiteOutput(capacity,T) for _=1:state.batch_size]
    else
        numerical=state.scratch.numerical; outputs=state.scratch.outputs
    end
    NI.update_prepared!(state.prepared,parameters)
    result=ChunkACEScratch(state.prepared,parameters,numerical,outputs,NL.normalize_properties(properties),zeros(Int,state.batch_size),zeros(Int,state.batch_size))
    state.scratch=result
    result
end
function create_scratch(c::ACECalculator;kwargs...)
    state=default_evaluation_state!(c)
    create_scratch(c,state;kwargs...)
end
# This extends the existing PFW evaluate! generic, preserving pair dispatch.
function evaluate!(buffers::NL.ChunkBuffers,s::ChunkACEScratch,c::ACECalculator,site)
    capacity=NI.neighbour_count(site)
    if length(s.numerical.Rs)<capacity
        T=eltype(s.prepared.weights)
        s.numerical=NI.make_scratch(s.prepared,capacity,T)
        for out in s.outputs; resize!(out.gradients,capacity); end
    end
    NI.evaluate_site!(s.outputs[1],s.prepared,site,s.parameters,s.numerical,s.properties)
    NL.scatter_site_output!(buffers,site,s.outputs[1],s.properties)
    nothing
end
function default_evaluation_state!(c::ACECalculator)
    isnothing(c.evaluation_state) && (c.evaluation_state=ACEEvaluationState(c))
    c.evaluation_state
end
function chunk_calculator_result(sys,c,state::ACEEvaluationState,parameters;
                                 energy_only=false,simd=true,fused=false,product_algorithm=state.prepared.product_algorithm,
                                 centers=nothing,properties=(AC.Energy(),AC.Forces(),AC.Virial()),parallel=false,kwargs...)
    state.prepared.model===c.model || throw(ArgumentError("evaluation state belongs to a different model"))
    state.active && throw(ArgumentError("ACE evaluation state already has an active call"))
    state.active=true
    try
        T=eltype(state.prepared.weights)
        fast=parameter_seed(parameters) isa Union{Float32,Float64}
        boundary=fast && T===Float32 ? precision_geometry(sys,Float32) : sys
        input=ace_inputs(boundary,c.model;environment=state.environment,precision=T)
        if !fast
            y=ace_evaluate(input.positions,c.model,parameters,input.environment,input.cell;state=state.model_state,energy_only,centers)
        else
            fused && !simd && throw(ArgumentError("tiled fusion requires SIMD radial blocks"))
            p=state.prepared
            state.prepared=typeof(p)(p.model,p.basis,p.weights,simd,fused,check_product_algorithm(product_algorithm),p.product_plan)
            n=NL.chunk_neighbours(boundary,state.environment,input.positions,input.cell;precision=T)
            state.chunks=NL.form_numa_workers(n,1;scatter_capacity=state.scatter_capacity)
            props=energy_only ? (AC.Energy(),) : NL.normalize_properties(properties)
            selected=falses(length(sys))
            selected[NL.selected_centers(centers,length(sys))].=true
            scratch=create_scratch(c,state;parameters,properties=props)
            scatter=NL.wants(AC.Forces,props) ? (forces=SVector{3,T},) : NamedTuple()
            sums=merge(NL.wants(AC.Energy,props) ? (energy=T,) : NamedTuple(),
                       NL.wants(AC.Virial,props) ? (virial=SMatrix{3,3,T,9},) : NamedTuple())
            run=() -> NL.foreach_worker!(state.chunks,1;scatter,sum_outputs=sums,init_scratch=()->scratch) do worker
                s=NL.get_scratch(worker)
                NL.evaluate_chunk_sites!(worker,s.prepared,s.parameters,s.numerical,s.outputs,s.ordinals,s.centers,s.properties,selected)
            end
            has_source(c.model) ? Base.invokelatest(run) : run()
            r=NL.get_results(state.chunks)
            forces=haskey(r,:forces) ? r.forces : zeros(SVector{3,T},length(sys))
            y=(energy=get(r,:energy,zero(T)),forces=reshape(collect(Iterators.flatten(forces)),3,:),virial=get(r,:virial,zero(SMatrix{3,3,T,9})))
        end
        energy_only && return (;energy=y.energy*u"eV")
        (;energy=y.energy*u"eV",forces=force_vectors(y.forces*u"eV/Å"),virial=static_virial(y.virial*u"eV"))
    finally
        state.active=false
    end
end

"""Prepare complete chunk neighborhoods with the calculator's explicit element
channel order (provider-native species IDs must not be used as ACE channels)."""
function NL.prepare_neighbors(sys,c::ACECalculator;precision=eltype(c.executor.force),width=c.environment.width,kwargs...)
    channels=[element_channel(c.model,AtomsBase.species(sys,i)) for i=1:length(sys)]
    NL.prepare_neighbors(sys,NI.search_cutoff(c.model);precision,width,channels,species_count=length(c.model.atomic_numbers),kwargs...)
end
