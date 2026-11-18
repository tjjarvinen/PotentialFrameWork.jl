module PotentialFrameWorkLuxCoreExt
using PotentialFrameWork, LuxCore
import AtomsCalculators as AC
import PotentialFrameWork.PairPotentials: lux_layer
using PotentialFrameWork.ACE: ResolvedModel, ace_evaluate, trainable_parameters

struct ChunkLayer{M,P} <: LuxCore.AbstractLuxLayer
    model::M
    initial::P
end
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
function lux_layer(m::PotentialFrameWork.ACE.ACEFeatureLayer;parameters=PotentialFrameWork.ACE.feature_parameters(m),executor=:functional)
    executor in (:functional,:chunks) || throw(ArgumentError("unknown feature executor"))
    executor===:chunks ? ChunkLayer(m,deepcopy(parameters)) : FeatureLayer(m,deepcopy(parameters))
end
function lux_layer(m::ResolvedModel;parameters,executor=:functional)
    executor in (:functional,:chunks) || throw(ArgumentError("unknown ACE executor"))
    executor===:chunks ? ChunkLayer(m,trainable_parameters(m,parameters)) : ACELayer(m,trainable_parameters(m,parameters))
end
LuxCore.initialparameters(rng::LuxCore.Random.AbstractRNG,l::Union{PairLayer,ACELayer,FeatureLayer,ChunkLayer})=deepcopy(l.initial)
LuxCore.initialstates(rng::LuxCore.Random.AbstractRNG,l::Union{PairLayer,ACELayer,FeatureLayer,ChunkLayer})=NamedTuple()
LuxCore.initialparameters(::Nothing,l::Union{PairLayer,ACELayer,FeatureLayer,ChunkLayer})=deepcopy(l.initial)
LuxCore.initialstates(::Nothing,l::Union{PairLayer,ACELayer,FeatureLayer,ChunkLayer})=NamedTuple()
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

function (l::ChunkLayer)(input,ps,st)
    A=PotentialFrameWork.ACE
    y=if l.model isa A.ACEFeatureLayer
        A.ace_features_chunks(input.positions,l.model,ps,input.neighbors;centers=get(input,:centers,nothing),batch_size=get(input,:batch_size,1),scatter_capacity=get(input,:scatter_capacity,4096))
    elseif l.model isa A.ACEFeatureEnergy
        A.feature_energy_evaluate(input.positions,l.model,ps,input.neighbors;centers=get(input,:centers,nothing),batch_size=get(input,:batch_size,1),scatter_capacity=get(input,:scatter_capacity,4096))
    else
        A.ace_evaluate_chunks(input.positions,l.model,ps,input.neighbors;centers=get(input,:centers,nothing),batch_size=get(input,:batch_size,1),scatter_capacity=get(input,:scatter_capacity,4096))
    end
    y,st
end

lux_layer(m::PotentialFrameWork.ACE.ACEFeatureEnergy;parameters)=ChunkLayer(m,deepcopy(parameters))

end
