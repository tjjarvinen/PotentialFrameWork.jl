module PotentialFrameWorkChainRulesCoreExt
using PotentialFrameWork, ChainRulesCore, Unitful
import AtomsCalculators as AC
import ChainRulesCore: rrule
using PotentialFrameWork.ACE: ace_evaluate, ace_vjp, ResolvedModel, FullEnvironment

function rrule(::typeof(parameter_energy),sys,p::AbstractPairPotential,parameters::NamedTuple;kwargs...)
    y=parameter_energy(sys,p,parameters;kwargs...)
    # Capture derivatives now: a supplied mutable workspace, system or parameter
    # vector may be reused/changed before the pullback is called.
    g=energy_parameter_gradient(sys,p,parameters;kwargs...)
    project=ProjectTo(parameters)
    function energy_pullback(delta)
        delta=unthunk(delta)
        dp=if delta isa AbstractZero
            ZeroTangent()
        else
            project(PotentialFrameWork.PairPotentials.parameter_scale(delta,g))
        end
        (NoTangent(),@not_implemented("system derivatives are outside parameter-only energy AD"),
         @not_implemented("differentiate explicit parameters, not the bundled potential"),dp)
    end
    y,energy_pullback
end

function rrule(::typeof(pair_evaluate),X,p::AbstractPairPotential,ps::NamedTuple,
               env::PairEnvironment,H;request=EnergyForcesVirial(),state=NamedTuple())
    y=pair_evaluate(X,p,ps,env,H;request,state)
    # Own every numerical input and discrete map needed for recomputation. Neither
    # caller mutation nor reuse of an unrelated inference workspace changes this.
    saved=deepcopy((X,p,ps,env,H))
    projectX=ProjectTo(X); projectps=ProjectTo(ps)
    function pair_pullback(delta)
        delta=unthunk(delta)
        if delta isa AbstractZero
            return (NoTangent(),ZeroTangent(),NoTangent(),ZeroTangent(),NoTangent(),NoTangent())
        end
        cot(name,zero_value)=if !hasproperty(delta,name)
            zero_value
        else
            v=unthunk(getproperty(delta,name))
            v isa AbstractZero ? zero_value : v
        end
        de=cot(:energy,zero(y.energy))
        df=cot(:forces,zeros(eltype(X),size(X)))
        dw=cot(:virial,zeros(eltype(X),3,3))
        g=pair_vjp(saved...,de,df,dw)
        (NoTangent(),projectX(g.positions),NoTangent(),projectps(g.parameters),NoTangent(),
         @not_implemented("cell derivatives are outside the fixed-cell pair AD contract"))
    end
    y,pair_pullback
end

function rrule(::typeof(ace_evaluate),X,model::ResolvedModel,ps,env::FullEnvironment,H;
               state=NamedTuple(),energy_only=false,centers=nothing)
    y=ace_evaluate(X,model,ps,env,H;state,energy_only,centers)
    saved_centers=deepcopy(centers)
    saved=deepcopy((X,model,ps,env,H))
    px=ProjectTo(X); pp=ProjectTo(ps)
    function ace_pullback(delta)
        delta=unthunk(delta)
        if delta isa AbstractZero
            return (NoTangent(),ZeroTangent(),NoTangent(),ZeroTangent(),NoTangent(),NoTangent())
        end
        function cot(name,z)
            hasproperty(delta,name) || return z
            v=unthunk(getproperty(delta,name))
            v isa AbstractZero ? z : v
        end
        de=cot(:energy,zero(y.energy)); df=cot(:forces,zeros(eltype(X),size(X))); dw=cot(:virial,zeros(eltype(X),3,3))
        g=ace_vjp(saved...,de,df,dw;centers=saved_centers)
        (NoTangent(),px(g.positions),NoTangent(),pp(g.parameters),NoTangent(),
         @not_implemented("cell derivatives are outside the fixed-cell ACE AD contract"))
    end
    y,ace_pullback
end
function rrule(::typeof(PotentialFrameWork.ACE.ace_features),X,
               layer::PotentialFrameWork.ACE.ACEFeatureLayer,ps,env::FullEnvironment,H)
    A=PotentialFrameWork.ACE
    y=A.ace_features(X,layer,ps,env,H)
    saved=deepcopy((X,layer,ps,env,H))
    px=ProjectTo(X); pp=ProjectTo(ps)
    function feature_pullback(delta)
        delta=unthunk(delta)
        if delta isa AbstractZero
            return (NoTangent(),ZeroTangent(),NoTangent(),ZeroTangent(),NoTangent(),NoTangent())
        end
        cot=map(eachindex(y)) do b
            d=unthunk(delta[b])
            d isa AbstractZero ? zero(y[b]) : d
        end
        all(size(cot[b])==size(y[b]) for b in eachindex(y)) || throw(DimensionMismatch("feature cotangent block shape"))
        g=A.feature_vjp(saved...,Tuple(cot))
        (NoTangent(),px(g.positions),NoTangent(),pp(g.parameters),NoTangent(),
         @not_implemented("cell derivatives are outside the fixed-cell feature AD contract"))
    end
    y,feature_pullback
end

# These entrypoints deliberately implement inference only.
for f in (PotentialFrameWork.PairPotentials.evaluate_gpu,
          PotentialFrameWork.PairPotentials.evaluate_gpu!,
          PotentialFrameWork.PairPotentials.execute_gpu!,
          PotentialFrameWork.PairPotentials.prepare_gpu_execution,
          PotentialFrameWork.NeighbourList.gpu_execute!)
    @eval function ChainRulesCore.rrule(::typeof($f),args...;kwargs...)
        throw(ArgumentError("GPU differentiation is unsupported; use the functional CPU pair_evaluate path"))
    end
    @eval function ChainRulesCore.frule(tangents,::typeof($f),args...;kwargs...)
        throw(ArgumentError("GPU differentiation is unsupported; use the functional CPU pair_evaluate path"))
    end
end


function rrule(::typeof(PotentialFrameWork.NeighbourList.site_evaluate),X,p,ps,n;input_derivatives=PotentialFrameWork.NeighbourInterface.SiteInputs(),kwargs...)
    NL=PotentialFrameWork.NeighbourList
    y=NL.site_evaluate(X,p,ps,n;kwargs...)
    # Save mathematical inputs only: no workgroup, worker scratch or caches.
    saved=deepcopy((X,p,ps,n)); options=deepcopy((;kwargs...))
    px=ProjectTo(X); pp=ProjectTo(ps)
    feature_templates=isnothing(y.features) ? nothing : [map(b -> (size(b),eltype(b)),h) for h in y.features]
    function pullback(delta)
        d=unthunk(delta)
        if d isa AbstractZero
            return (NoTangent(),ZeroTangent(),NoTangent(),ZeroTangent(),NoTangent())
        end
        function normalize(v,templates)
            v=unthunk(v)
            (v isa AbstractZero || isnothing(v)) && return nothing
            map(eachindex(templates)) do k
                b=unthunk(v[k]); shape,T=templates[k]
                (b isa AbstractZero || isnothing(b)) ? zeros(T,shape) : b
            end |> Tuple
        end
        function cot(key)
            hasproperty(d,key) || return nothing
            v=unthunk(getproperty(d,key))
            v isa AbstractZero && return nothing
            if key===:features
                if isnothing(feature_templates)
                    NL.nonzero_bar(v) && throw(ArgumentError("feature cotangent supplied for absent outputs"))
                    return nothing
                end
                return [normalize(v[i],feature_templates[i]) for i in eachindex(feature_templates)]
            end
            v
        end
        bars=(energy=cot(:energy),forces=cot(:forces),virial=cot(:virial),features=cot(:features))
        g=NL.site_pullback(saved...,bars;request=input_derivatives,options...)
        (NoTangent(),isnothing(g.positions) ? ZeroTangent() : px(g.positions),NoTangent(),isnothing(g.parameters) ? ZeroTangent() : pp(g.parameters),NoTangent())
    end
    y,pullback
end

function rrule(::typeof(PotentialFrameWork.ACE.ace_features_chunks),X,l,ps,n;kwargs...)
    A=PotentialFrameWork.ACE; NL=PotentialFrameWork.NeighbourList; NI=PotentialFrameWork.NeighbourInterface
    r,back=rrule(NL.site_evaluate,X,A.FeatureSiteModel(l),ps,n;outputs=NI.SiteOutputs(false,false,true),kwargs...)
    N=size(X,2)
    y=A.stack_feature_sites(r.features,l,N)
    shapes=map(size,y); types=map(eltype,y)
    function pullback(delta)
        d=unthunk(delta)
        d isa AbstractZero && return back(ZeroTangent())
        blocks=map(eachindex(y)) do b
            v=unthunk(d[b]); v isa AbstractZero ? zeros(types[b],shapes[b]) : v
        end
        bars=[Tuple(copy(v[:,:,i]) for v in blocks) for i in 1:N]
        back((features=bars,))
    end
    y,pullback
end
function rrule(::typeof(PotentialFrameWork.ACE.ace_evaluate_chunks),X,m,ps,n;kwargs...)
    A=PotentialFrameWork.ACE; NL=PotentialFrameWork.NeighbourList
    r,back=rrule(NL.site_evaluate,X,A.ScalarSiteModel(m),ps,n;kwargs...)
    (;energy=r.energy,forces=r.forces,virial=r.virial),back
end

# Unitful's numerical magnitude is the trainable scalar. Zygote represents
# ustrip cotangents structurally as Tangent(val=...), including SArray elements.
function unit_scalar_bar(d,unit)
    d=unthunk(d)
    (d isa AbstractZero || isnothing(d)) && return 0.0
    d isa Unitful.AbstractQuantity && return ustrip(Unitful.NoUnits,d*unit)
    hasproperty(d,:val) && return unit_scalar_bar(unthunk(getproperty(d,:val)),unit)
    d isa Real && return d
    throw(ArgumentError("unsupported unitful scalar cotangent $(typeof(d))"))
end
function unit_component_bar(d,k)
    d=unthunk(d)
    (d isa AbstractZero || isnothing(d)) && return ZeroTangent()
    d isa Tangent && hasproperty(d,:data) && return unit_component_bar(getproperty(d,:data),k)
    unthunk(d[k])
end

# Explicit-parameter calculator training bypasses mutable inference state in its
# derivative rule. Unitful structural cotangents are normalized at this boundary.
function rrule(config::RuleConfig{>:HasReverseMode},::typeof(AC.calculate),props,sys,
               c::PotentialFrameWork.ACE.ACECalculator,ps,state;kwargs...)
    A=PotentialFrameWork.ACE; NL=PotentialFrameWork.NeighbourList; NI=PotentialFrameWork.NeighbourInterface
    normalized=NL.normalize_properties(props)
    if state isa A.ACEEvaluationState
        state.active && throw(ArgumentError("ACE evaluation state already has an active call"))
        state.prepared.model===c.model || throw(ArgumentError("evaluation state belongs to a different model"))
    end
    # Preserve inference keyword validation, state updates, precision and public
    # units exactly. The separately owned numerical rule supplies the backward.
    y=AC.calculate(props,sys,c,ps,state;kwargs...)
    state isa A.ACEEvaluationState && (state.active=true)
    r,back=try
        n=NL.prepare_neighbors(sys,c)
        options=(;kwargs...)
        centers=get(options,:centers,nothing)
        outputs=NI.SiteOutputs(NL.wants(AC.Energy,normalized),NL.wants(AC.Forces,normalized)||NL.wants(AC.Virial,normalized),false)
        rrule(NL.site_evaluate,n.positions,A.ScalarSiteModel(c.model),ps,n;outputs,centers,
            input_derivatives=NI.SiteInputs(false,true),batch_size=c.executor.batch_size,scatter_capacity=state isa A.ACEEvaluationState ? state.scatter_capacity : 4096)
    finally
        state isa A.ACEEvaluationState && (state.active=false)
    end
    N=length(sys)
    function pullback(d)
        d=unthunk(d)
        d isa AbstractZero && return (NoTangent(),NoTangent(),NoTangent(),NoTangent(),ZeroTangent(),NoTangent())
        field(k)=hasproperty(d,k) ? unthunk(getproperty(d,k)) : ZeroTangent()
        e=field(:energy); f=field(:forces); w=field(:virial)
        numeric=(energy=unit_scalar_bar(e,u"eV"),
            forces=[unit_scalar_bar(unit_component_bar(unit_component_bar(f,i),a),u"eV/Å") for a=1:3,i=1:N],
            virial=[unit_scalar_bar(unit_component_bar(w,k),u"eV") for k=1:9] |> x -> reshape(x,3,3))
        g=back(numeric)
        (NoTangent(),NoTangent(),NoTangent(),NoTangent(),g[4],NoTangent())
    end
    y,pullback
end
function rrule(config::RuleConfig{>:HasReverseMode},::typeof(AC.calculate),props,sys,
               c::PotentialFrameWork.ACE.ACECalculator,ps;state=c.state,kwargs...)
    y,back=rrule(config,AC.calculate,props,sys,c,ps,state;kwargs...)
    y,d -> back(d)[1:5]
end

function rrule(::typeof(PotentialFrameWork.ACE.feature_energy_evaluate),X,m,ps,n;kwargs...)
    A=PotentialFrameWork.ACE; NL=PotentialFrameWork.NeighbourList
    rrule(NL.site_evaluate,X,A.FeatureEnergySiteModel(m),ps,n;kwargs...)
end

end
