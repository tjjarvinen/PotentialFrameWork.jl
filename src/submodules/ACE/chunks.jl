"""Exclusive mutable CPU inference state for serial or threaded site chunks.

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
    worker_scratch::Vector{Any}
    layout_valid::Bool
    periodic::NTuple{3,Bool}
    scatter_capacity::Int
    batch_size::Int
    active::Bool
    guard::ReentrantLock
end
function ACEEvaluationState(c::ACECalculator;state=deepcopy(c.state),scatter_capacity=4096)
    scatter_capacity isa Integer && scatter_capacity>0 || throw(ArgumentError("scatter_capacity must be positive"))
    # The angular provider, derived weights and product workspace are state-owned.
    original=c.executor.prepared
    p=typeof(original)(c.model,deepcopy(original.basis),copy(original.weights),original.simd,original.fused,original.product_algorithm,original.product_plan)
    ACEEvaluationState(FullEnvironment(;width=c.environment.width),p,state,nothing,nothing,Any[],false,(false,false,false),Int(scatter_capacity),c.executor.batch_size,false,ReentrantLock())
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
    NI.update_prepared!(state.prepared,parameters)
    state.scratch=chunk_ace_scratch(state.prepared,parameters,NL.normalize_properties(properties),
                                    state.scratch,chunk_site_capacity(state),state.batch_size)
end

chunk_site_capacity(state)=maximum((state.environment.offsets[end,i]-state.environment.offsets[1,i]
                                    for i in eachindex(state.environment.species));init=0)
function chunk_ace_scratch(prepared,ps,props,old,capacity,batch_size)
    T=eltype(prepared.weights)
    if isnothing(old) || length(old.numerical.Rs)<capacity
        numerical=NI.make_scratch(prepared,capacity,T)
        outputs=[NL.SiteOutput(capacity,T) for _ in 1:batch_size]
    else
        numerical=old.numerical
        outputs=old.outputs
    end
    ordinals=isnothing(old) ? zeros(Int,batch_size) : old.ordinals
    centers=isnothing(old) ? zeros(Int,batch_size) : old.centers
    ChunkACEScratch(prepared,ps,numerical,outputs,props,ordinals,centers)
end

function worker_ace_scratch!(state,slot,ps,props,capacity)
    old=state.worker_scratch[slot]
    p=state.prepared
    basis=isnothing(old) ? deepcopy(p.basis) : old.prepared.basis
    prepared=typeof(p)(p.model,basis,p.weights,p.simd,p.fused,p.product_algorithm,p.product_plan)
    state.worker_scratch[slot]=chunk_ace_scratch(prepared,ps,props,old,capacity,state.batch_size)
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
const ACE_STATE_INIT_LOCK=ReentrantLock()
function default_evaluation_state!(c::ACECalculator)
    lock(ACE_STATE_INIT_LOCK) do
        isnothing(c.evaluation_state) && (c.evaluation_state=ACEEvaluationState(c))
        c.evaluation_state
    end
end

function acquire_ace_state!(state)
    NL.acquire_chunk_access!(state.guard,()->state.active,"ACE evaluation state already has an active call")
end

"""Join the state's idle workers and release their scratch. The state is reusable."""
function Base.close(state::ACEEvaluationState)
    acquire_ace_state!(state)
    try
        isnothing(state.chunks) || close(state.chunks)
        empty!(state.worker_scratch)
    finally
        unlock(state.guard)
    end
    nothing
end

function chunk_geometry!(state,sys,model,T,groups,reuse_topology)
    old=state.chunks
    periodic=Tuple(Bool.(AtomsBase.periodicity(sys)))
    if reuse_topology
        state.layout_valid || throw(ArgumentError("prepare a fresh chunk layout before reusing topology"))
        min(groups,length(old.neighbours.ranges))==length(old.groups) &&
            state.scatter_capacity==old.scatter_capacity || throw(ArgumentError("chunk layout options changed; rebuild topology"))
        X,H=ace_geometry(sys)
        periodic==state.periodic && H==old.neighbours.cell ||
            throw(ArgumentError("cell or periodicity changed; rebuild topology"))
        channels=[element_channel(model,AtomsBase.species(sys,i)) for i in 1:length(sys)]
        channels==state.environment.species || throw(ArgumentError("atom count or species changed; rebuild topology"))
        NL.refresh_positions!(old,X)
        return old
    end
    state.layout_valid=false
    if !isnothing(old)
        empty!(old.fields)
        old.results=nothing
    end
    input=ace_inputs(sys,model;environment=state.environment,precision=T)
    n=NL.chunk_neighbours(sys,state.environment,input.positions,input.cell;precision=T)
    fresh=NL.form_numa_workers(n,groups;scatter_capacity=state.scatter_capacity)
    if !isnothing(old)
        fresh.pool=old.pool
        old.pool=nothing
    end
    state.chunks=fresh
    state.periodic=periodic
    state.layout_valid=true
    fresh
end

function run_ace_chunks!(w,scratch,selected,workers,scatter,sums,timings)
    NL.foreach_worker!(w,workers;scatter,sum_outputs=sums,slot_scratch=true,
                       init_scratch=slot -> scratch[slot],timings) do worker
        s=NL.get_scratch(worker)
        NL.evaluate_chunk_sites!(worker,s.prepared,s.parameters,s.numerical,s.outputs,s.ordinals,s.centers,s.properties,selected)
    end
end

function chunk_calculator_result(sys,c,state::ACEEvaluationState,parameters;
                                 energy_only=false,simd=true,fused=false,product_algorithm=state.prepared.product_algorithm,
                                 centers=nothing,properties=(AC.Energy(),AC.Forces(),AC.Virial()),parallel=false,
                                 workers=1,groups=1,reuse_topology=false,timings=nothing,kwargs...)
    state.prepared.model===c.model || throw(ArgumentError("evaluation state belongs to a different model"))
    acquire_ace_state!(state)
    state.active=true
    try
        T=eltype(state.prepared.weights)
        fast=parameter_seed(parameters) isa Union{Float32,Float64}
        boundary=fast && T===Float32 ? precision_geometry(sys,Float32) : sys
        if !fast
            reuse_topology && throw(ArgumentError("topology reuse requires numerical ACE parameters"))
            state.layout_valid=false
            input=ace_inputs(boundary,c.model;environment=state.environment,precision=T)
            y=ace_evaluate(input.positions,c.model,parameters,input.environment,input.cell;state=state.model_state,energy_only,centers)
        else
            fused && !simd && throw(ArgumentError("tiled fusion requires SIMD radial blocks"))
            p=state.prepared
            state.prepared=typeof(p)(p.model,p.basis,p.weights,simd,fused,check_product_algorithm(product_algorithm),p.product_plan)
            w=chunk_geometry!(state,boundary,c.model,T,groups,reuse_topology)
            props=energy_only ? (AC.Energy(),) : NL.normalize_properties(properties)
            selected=falses(length(sys))
            selected[NL.selected_centers(centers,length(sys))].=true
            scatter=NL.wants(AC.Forces,props) ? (forces=SVector{3,T},) : NamedTuple()
            sums=merge(NL.wants(AC.Energy,props) ? (energy=T,) : NamedTuple(),
                       NL.wants(AC.Virial,props) ? (virial=SMatrix{3,3,T,9},) : NamedTuple())
            scratch=if workers==1 || length(w.chunks)==1
                [create_scratch(c,state;parameters,properties=props)]
            else
                NI.update_prepared!(state.prepared,parameters)
                capacity=chunk_site_capacity(state)
                slots=min(workers,length(w.chunks))
                while length(state.worker_scratch)<slots; push!(state.worker_scratch,nothing); end
                if slots==0
                    ()
                else
                    first_scratch=worker_ace_scratch!(state,1,parameters,props,capacity)
                    entries=Vector{typeof(first_scratch)}(undef,slots)
                    entries[1]=first_scratch
                    for slot in 2:slots
                        entries[slot]=worker_ace_scratch!(state,slot,parameters,props,capacity)
                    end
                    entries
                end
            end
            if has_source(c.model)
                Base.invokelatest(run_ace_chunks!,w,scratch,selected,workers,scatter,sums,timings)
            else
                run_ace_chunks!(w,scratch,selected,workers,scatter,sums,timings)
            end
            export_start=NL.chunk_clock(timings)
            r=NL.get_results(w)
            energy=get(r,:energy,zero(T))*u"eV"
            result=if energy_only
                (;energy)
            else
                raw_forces=get(() -> zeros(SVector{3,T},length(sys)),r,:forces)
                forces=raw_forces .* u"eV/Å"
                virial=get(r,:virial,zero(SMatrix{3,3,T,9}))*u"eV"
                (;energy,forces,virial)
            end
            NL.record_chunk_export!(timings,export_start)
            return result
        end
        energy_only && return (;energy=y.energy*u"eV")
        (;energy=y.energy*u"eV",forces=force_vectors(y.forces*u"eV/Å"),virial=static_virial(y.virial*u"eV"))
    finally
        state.active=false
        unlock(state.guard)
    end
end

"""Prepare complete chunk neighborhoods with the calculator's explicit element
channel order (provider-native species IDs must not be used as ACE channels)."""
function NL.prepare_neighbors(sys,c::ACECalculator;precision=eltype(c.executor.force),width=c.environment.width,kwargs...)
    channels=[element_channel(c.model,AtomsBase.species(sys,i)) for i=1:length(sys)]
    NL.prepare_neighbors(sys,NI.search_cutoff(c.model);precision,width,channels,species_count=length(c.model.atomic_numbers),kwargs...)
end
