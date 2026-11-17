# Run with: julia --project examples/site_chunks.jl
using PotentialFrameWork, AtomsBase, Unitful, StaticArrays
using PotentialFrameWork.NeighbourList
const NI=PotentialFrameWork.NeighbourInterface

system=AtomsBase.atomic_system([:H=>[3.1,0.,0.]*u"Å",:He=>[0.2,0.,0.]*u"Å",
    :H=>[1.4,0.3,0.]*u"Å"],ntuple(a -> SVector{3,Float64}(ntuple(b -> a==b ? 8. : 0.,3))*u"Å",3),(true,false,false))
workgroups=form_numa_workers(prepare_neighbors(system,2.4),1;scatter_capacity=3)

# First pass writes one model-defined feature per center, directly in cell order.
foreach_worker!(workgroups,1;per_site=(layer1=SVector{2,Float64},)) do worker
    buffers=get_buffers(worker)
    for site in worker
        i=NI.central_index(site) # original index, independent of cell ordering
        buffers.per_site.layer1[site_ordinal(site)]=SVector(Float64(i),0.1i^2)
    end
end

# A second pass reads both center and neighboring features. Multiple periodic
# images read the same atom feature, retaining their individual geometry.
foreach_worker!(workgroups,1;site_inputs=(:layer1,),per_site=(layer2=SVector{2,Float64},)) do worker
    buffers=get_buffers(worker)
    for site in worker
        feature=central_feature(site,Val(:layer1))
        for edge=1:NI.neighbour_count(site)
            NI.active(site,edge) || continue
            d=NI.displacement(site,edge)
            feature+=neighbor_feature(site,edge,Val(:layer1))*exp(-sum(abs2,d))
        end
        buffers.per_site.layer2[site_ordinal(site)]=feature
    end
end
println("layer2 in original atom order: ",export_field(workgroups,:layer2))
