module PotentialFrameWorkLuxCoreExt
using PotentialFrameWork, LuxCore
import AtomsCalculators as AC
import PotentialFrameWork.PairPotentials: lux_layer
using PotentialFrameWork.ACE: ResolvedModel, ace_evaluate, trainable_parameters

struct PairLayer{M,P} <: LuxCore.AbstractLuxLayer
    model::M
    initial::P
end
struct FeatureLayer{M,P} <: LuxCore.AbstractLuxLayer
    model::M
    initial::P
end
struct ACELayer{M,P} <: LuxCore.AbstractLuxLayer
    model::M
    initial::P
end
lux_layer(p::AbstractPairPotential;parameters=AC.get_parameters(p))=PairLayer(p,deepcopy(parameters))
lux_layer(m::PotentialFrameWork.ACE.ACEFeatureLayer;parameters=PotentialFrameWork.ACE.feature_parameters(m))=FeatureLayer(m,deepcopy(parameters))
lux_layer(m::ResolvedModel;parameters)=ACELayer(m,trainable_parameters(m,parameters))
LuxCore.initialparameters(rng::LuxCore.Random.AbstractRNG,l::Union{PairLayer,ACELayer,FeatureLayer})=deepcopy(l.initial)
LuxCore.initialstates(rng::LuxCore.Random.AbstractRNG,l::Union{PairLayer,ACELayer,FeatureLayer})=NamedTuple()
LuxCore.initialparameters(::Nothing,l::Union{PairLayer,ACELayer,FeatureLayer})=deepcopy(l.initial)
LuxCore.initialstates(::Nothing,l::Union{PairLayer,ACELayer,FeatureLayer})=NamedTuple()
function (l::PairLayer)(input,ps,st)
    y=pair_evaluate(input.positions,l.model,ps,input.environment,input.cell)
    y,st
end
function (l::ACELayer)(input,ps,st)
    y=ace_evaluate(input.positions,l.model,ps,input.environment,input.cell;centers=get(input,:centers,nothing))
    y,st
end
function (l::FeatureLayer)(input,ps,st)
    y=PotentialFrameWork.ACE.ace_features(input.positions,l.model,ps,input.environment,input.cell)
    y,st
end

end
