# julia --project=<test/application-env> --threads=4 benchmark/neighbour_evaluation.jl <package-root> <report.toml>
# Run before/after sequentially with identical inputs and no concurrent tests.
# package-root can be a preserved baseline checkout; no source/model is modified.
using PotentialFrameWork, AtomsBase, Unitful, SpheriCart, LinearAlgebra, TOML, SIMD
@assert startswith(realpath(pathof(PotentialFrameWork)),realpath(ARGS[1])) "activate an environment pointing to the requested source"
const A=PotentialFrameWork.ACE
const NL=PotentialFrameWork.NeighbourList
const PP=PotentialFrameWork.PairPotentials
BLAS.set_num_threads(1)
model,ps,st=A.load_resolved(joinpath(ARGS[1],"data/ace_ref/resolved_model_v1.toml"))
sys=AtomsBase.atomic_system([:Ti=>[0.1,0.2,0.3]*u"Å",:Al=>[2.3,0.7,1.2]*u"Å",:Ti=>[1.1,2.4,1.5]*u"Å"],
    ([4.7,0.,0.]*u"Å",[0.5,4.8,0.]*u"Å",[0.1,0.3,4.6]*u"Å"),(true,true,true))
input=A.ace_inputs(sys,model;environment=A.FullEnvironment(width=4))
w=A.CPUExecutor(model,input.environment;slots=4)
p=SmoothPair([:Ti,:Al],[1.0 0.7;0.7 0.4],3.2)
pw=Workspace(width=4); NL.prepare!(pw,sys,cutoff(p)); bound=bind_species(p,pw.species_map)
function measure(f)
    for _=1:20; f(); end
    times=[(@elapsed begin; for _=1:100; f(); end; end)/100 for s=1:9]
    (;seconds_min=minimum(times),seconds_median=sort(times)[5],bytes=@allocated(f()))
end
report=Dict{String,Any}("julia"=>string(VERSION),"cpu"=>Sys.cpu_info()[1].model,"threads"=>Threads.nthreads(),"blas_threads"=>BLAS.get_num_threads(),"SIMD"=>string(pkgversion(SIMD)),"SpheriCart"=>string(pkgversion(SpheriCart)),"AtomsBase"=>string(pkgversion(AtomsBase)),"Unitful"=>string(pkgversion(Unitful)),"method"=>"serial prepared E/F/W; width 4; four slots; BLAS 1; 20 warmups; nine samples of 100 calls; three periodic atoms; preparation/output excluded")
report["pair"]=Dict(pairs(measure(()->PP.execute!(pw,bound,EnergyForcesVirial()))))
if isdefined(PP,:PreparedSmooth)
    NI=PotentialFrameWork.NeighbourInterface
    native=NI.prepare_potential(p,pw.species_map)
    props=(NI.AtomsCalculators.Energy(),NI.AtomsCalculators.Forces(),NI.AtomsCalculators.Virial())
    report["native_radial_smooth"]=Dict(pairs(measure(()->NL.execute!(pw,native,props))))
    radial=(r,si,sj,req) -> req isa EnergyOnly ? (;energy=(1-r*r/3.2^2)^3) :
        (;energy=(1-r*r/3.2^2)^3,derivative=-6r/3.2^2*(1-r*r/3.2^2)^2)
    adapter=NI.prepare_potential(RadialPotential(radial,3.2),pw.species_map)
    report["native_radial_adapter"]=Dict(pairs(measure(()->NL.execute!(pw,adapter,props))))
end
report["ace"]=Dict(pairs(measure(()->A.execute_prepared!(w,input.positions,input.cell,ps))))
report["scratch_bytes"]=Base.summarysize(w.scratch)
if hasproperty(w,:slots)
    report["slot_bytes"]=Base.summarysize(w.slots)
    report["gradient_bytes"]=sum(sizeof(o.gradients) for slot in w.slots for o in slot.outputs)
    for batch in (2,5)
        wb=A.CPUExecutor(model,input.environment;slots=4,batch_size=batch)
        report["ace_batch_$batch"]=Dict(pairs(measure(()->A.execute_prepared!(wb,input.positions,input.cell,ps))))
    end
end
open(io->TOML.print(io,report),ARGS[2],"w")
println(report)
