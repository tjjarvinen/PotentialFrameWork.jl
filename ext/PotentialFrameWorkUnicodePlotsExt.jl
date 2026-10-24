module PotentialFrameWorkUnicodePlotsExt
using PotentialFrameWork, UnicodePlots
using PotentialFrameWork.ACE: ResolvedModel, basis_summary, radial_samples
import PotentialFrameWork.ACE: plot_basis

"""Terminal plots from owning ACE inspection arrays. Coefficient magnitudes are
normalization/representation dependent; these plots do not rank importance.
Kinds: :orders, :weights, :coefficients, :radial, :coupling, :counts.
"""
function plot_basis(m::ResolvedModel,ps=nothing;kind=:orders,indices=1:min(64,size(m.coupling,1)),
                    channels=1:min(6,length(m.radial_spec)),parameters=NamedTuple(),kwargs...)
    summary=basis_summary(m,ps)
    if kind in (:orders,:weights)
        counts=kind===:orders ? summary.order_counts : summary.index_weight_counts
        labels=sort!(collect(keys(counts)))
        isempty(labels) && return barplot(["empty/undefined"],[0];kwargs...)
        return barplot(string.(labels),[counts[k] for k in labels];
            xlabel=kind===:orders ? "invariants per correlation order" : "invariants per index weight",kwargs...)
    elseif kind===:coefficients
        isnothing(ps) && throw(ArgumentError("provide parameters for coefficient distribution"))
        values=abs.(vec(ps.WB[collect(indices),:]))
        isempty(values) && return barplot(["no coefficients"],[0];kwargs...)
        return histogram(values;xlabel="|coefficient| (representation dependent)",kwargs...)
    elseif kind===:coupling
        C=m.coupling[collect(indices),:]
        isempty(C) && return barplot(["no coupling entries"],[0];kwargs...)
        return spy(C;kwargs...)
    elseif kind===:radial
        isempty(channels) && return barplot(["no radial channels"],[0];kwargs...)
        sample=radial_samples(m;channels,parameters)
        p=lineplot(sample.distances,sample.values[:,1];name="channel $(first(channels))",xlabel="r / Å",kwargs...)
        for j=2:length(channels)
            lineplot!(p,sample.distances,sample.values[:,j];name="channel $(channels[j])")
        end
        return p
    elseif kind===:counts
        ds=summary.descriptors[collect(indices)]
        os=sort!(unique(vcat(Int[],[d.orders for d in ds]...)))
        ws=sort!(unique([d.index_weight for d in ds if !isnothing(d.index_weight)]))
        (isempty(os) || isempty(ws)) && return barplot(["no grouped counts"],[0];kwargs...)
        counts=[count(d -> o in d.orders && d.index_weight==w,ds) for o in os,w in ws]
        return heatmap(counts;xlabel="index-weight groups $(ws)",ylabel="order groups $(os)",kwargs...)
    end
    throw(ArgumentError("unknown ACE plot kind $kind"))
end
end
