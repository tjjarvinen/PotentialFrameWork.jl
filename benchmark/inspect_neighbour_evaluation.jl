# julia --project=<application-env> benchmark/inspect_neighbour_evaluation.jl [warntype-output.txt]
using PotentialFrameWork, AtomsBase, Unitful, SpheriCart, InteractiveUtils
import AtomsCalculators as AC
const A=PotentialFrameWork.ACE
const NL=PotentialFrameWork.NeighbourList
const NI=PotentialFrameWork.NeighbourInterface
model,ps,st=A.load_resolved(joinpath(@__DIR__,"../data/ace_ref/resolved_model_v1.toml"))
sys=AtomsBase.atomic_system([:Ti=>[0.1,0.2,0.3]*u"Å",:Al=>[2.3,0.7,1.2]*u"Å",:Ti=>[1.1,2.4,1.5]*u"Å"],
    ([4.7,0.,0.]*u"Å",[0.5,4.8,0.]*u"Å",[0.1,0.3,4.6]*u"Å"),(true,true,true))
input=A.ace_inputs(sys,model;environment=A.FullEnvironment(width=4))
w=A.CPUExecutor(model,input.environment;slots=4)
A.execute_prepared!(w,input.positions,input.cell,ps)
p=w.prepared; site=NL.SiteView(w.environment,input.positions,input.cell,1)
out=w.slots[1].outputs[1]; props=(AC.Energy(),AC.Forces(),AC.Virial())
args=(out,p,site,ps,w.scratch[1],props)
NI.evaluate_site!(args...)
println("callback return types: ",Base.return_types(NI.evaluate_site!,Tuple{map(typeof,args)...}))
open(isempty(ARGS) ? "/tmp/pfw-neighbour-site-warntype.txt" : first(ARGS),"w") do io
    code_warntype(io,NI.evaluate_site!,Tuple{map(typeof,args)...};debuginfo=:none)
end
function callback_allocation(args)
    NI.evaluate_site!(args...)
    @allocated NI.evaluate_site!(args...)
end
println("function-barrier site callback bytes: ",callback_allocation(args))
pair=SmoothPair([:Ti,:Al],[1.0 0.7;0.7 0.4],3.2)
pair_workspace=Workspace(width=4)
NL.prepare!(pair_workspace,sys,cutoff(pair))
bound=PotentialFrameWork.PairPotentials.LegacyPair(bind_species(pair,pair_workspace.species_map))
pair_args=(pair_workspace,bound,props,1)
function pair_allocation(args)
    NL.cell_task!(args...)
    @allocated NL.cell_task!(args...)
end
println("function-barrier pair cell-task bytes: ",pair_allocation(pair_args))
open(length(ARGS)<2 ? "/tmp/pfw-neighbour-pair-warntype.txt" : ARGS[2],"w") do io
    code_warntype(io,NL.cell_task!,Tuple{map(typeof,pair_args)...};debuginfo=:none)
end
native=NI.prepare_potential(pair,pair_workspace.species_map)
native_args=(pair_workspace,native,props,1)
println("function-barrier native radial pair cell-task bytes: ",pair_allocation(native_args))
open("/tmp/pfw-radial-pair-warntype.txt","w") do io
    code_warntype(io,NL.cell_task!,Tuple{map(typeof,native_args)...};debuginfo=:none)
end
