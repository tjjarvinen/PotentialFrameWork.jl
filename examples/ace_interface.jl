# Application environment: PotentialFrameWork, AtomsBase, SpheriCart, Unitful.
using PotentialFrameWork, AtomsBase, SpheriCart, Unitful
import AtomsCalculators as AC
const A=PotentialFrameWork.ACE
model,ps,st=A.load_resolved(joinpath(@__DIR__,"../data/ace_ref/resolved_model_v1.toml"))
println(A.basis_summary(model).dimensions)
selection=A.BasisSelection(model;orders=1:2,max_index_weight=8)
masked=A.mask_basis(model,selection) # reversible; ps retains all original values
reduced=A.restrict_basis(model,ps,selection) # no executor allocated before reduction
projected=A.project_species(reduced.model,[:Ti],1:2)
sys=AtomsBase.isolated_system([:Ti=>[0.,0.,0.]*u"Å",:Al=>[2.3,0.7,1.2]*u"Å",:Ti=>[1.1,2.4,1.5]*u"Å"])
calc=A.ACECalculator(projected,reduced.parameters;slots=1)
println(AC.energy_forces_virial(sys,calc;centers=[1,3]).energy)
input=A.ace_inputs(sys,projected)
design=A.design_matrix(input.positions,projected,reduced.parameters,input.environment,input.cell;centers=[1,3])
println("Design dimensions: ",size(design.matrix))

# A trusted, explicit Julia source bundle. The factory returns a generic
# callable with separate trainable ps; use a registration instead for a package
# implementation. Saved source loading always requires allow_code=true.
metadata=A.RadialMetadata(model.radial_spec;cutoff=5.5,index_weights=sum.(model.radial_spec),
    simd=true,ad=true,origin=:regular)
config=(powers=[n-1 for (n,l) in model.radial_spec],cutoff=5.5)
source=read(joinpath(@__DIR__,"ace_radial_source.jl"),String)
channel=A.source_radial(source,"radial_factory",metadata;config,allow_code=true)
custom=A.replace_radial(model,ps,fill(channel,2,2);parameters=(rate=[0.7],amplitude=[0.2]))
println(A.radial_capabilities(channel,custom.parameters.custom_radial))
customcalc=A.ACECalculator(custom.model,custom.parameters;slots=1)
mktempdir() do dir
    path=joinpath(dir,"custom.toml")
    A.save_calculator(path,customcalc)
    restored=A.ACECalculator(path;allow_code=true,slots=1)
    println(AC.potential_energy(sys,restored))
end
# Optional, when UnicodePlots is installed/loaded:
# using UnicodePlots
# display(A.plot_basis(projected,reduced.parameters;kind=:orders))
