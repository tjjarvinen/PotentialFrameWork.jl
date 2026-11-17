# Storage/allocation diagnostic, not a throughput comparison with coloring.
# julia --project=<test-env> --threads=4 benchmark/cpu_chunks/memory.jl
using PotentialFrameWork, AtomsBase, Unitful, StaticArrays, LinearAlgebra
using PotentialFrameWork.NeighbourList
using Pkg
include(joinpath(@__DIR__,"../../test/fixtures.jl"))
const NI=PotentialFrameWork.NeighbourInterface
function pass!(w,temporary_bytes)
    foreach_worker!(w,1;scatter=(forces=SVector{3,Float64},),sum_outputs=(energy=Float64,)) do worker
        b=get_buffers(worker)
        temporary_bytes[]=max(temporary_bytes[],sizeof(parent(b.scatter.forces)))
        for site in worker
            b.sum_outputs.energy[]+=1
            for edge=1:NI.neighbour_count(site)
                NI.active(site,edge) || continue
                g=NI.displacement(site,edge)
                b.scatter.forces[central_slot(site)]+=g
                b.scatter.forces[neighbor_slot(site,edge)]-=g
            end
        end
    end
end
function measure()
    sys=TestSystem([1.4 1.5 0.5 2.4;0 0 0 0;0 0 0 0],10.0Matrix(I,3,3),(true,false,false),fill(:H,4))
    w=form_numa_workers(prepare_neighbors(sys,1.;width=4),1;scatter_capacity=3)
    temporary_bytes=Ref(0)
    pass!(w,temporary_bytes); pass!(w,temporary_bytes)
    bytes=@allocated pass!(w,temporary_bytes)
    println("Julia=",VERSION," threads=",Threads.nthreads()," CPU=",Sys.CPU_NAME)
    println("atoms=",length(sys)," cells=",length(w.neighbours.ranges)," chunks=",length(w.chunks)," max_destinations=",w.observed_capacity," configured_capacity=",w.scatter_capacity)
    println("Float64 force final bytes=",sizeof(w.results.forces)," bounded temporary bytes=",temporary_bytes[]," warmed pass allocated bytes=",bytes)
    for name in ("AtomsBase","AtomsCalculators","StaticArrays","SIMD","SpheriCart","KernelAbstractions","AcceleratedKernels","ForwardDiff","ChainRulesCore","Zygote","LuxCore")
        packages=filter(x -> x.name==name,collect(values(Pkg.dependencies())))
        isempty(packages) || println(name,"=",only(packages).version)
    end
end
measure()
