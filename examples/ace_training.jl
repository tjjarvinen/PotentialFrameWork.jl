# Application environment: PotentialFrameWork, AtomsBase, AtomsCalculators,
# SpheriCart, ACEfit, LuxCore, ForwardDiff, Zygote and Unitful.
using PotentialFrameWork, AtomsBase, SpheriCart, ACEfit, LuxCore, ForwardDiff, Zygote, Unitful
import AtomsCalculators as AC
const A=PotentialFrameWork.ACE

model,readout,st=A.load_recipe(joinpath(@__DIR__,"../data/ace_ref/reference_model.json"))
sys=AtomsBase.atomic_system([
    :Ti=>[0.1,0.2,0.3]*u"Å", :Al=>[2.5,0.7,1.2]*u"Å", :Ti=>[1.1,2.4,1.5]*u"Å"],
    ([9.,0.,0.]*u"Å",[0.,9.,0.]*u"Å",[0.,0.,9.]*u"Å"),(false,false,false))
input=A.ace_inputs(sys,model)
calc=A.ACECalculator(model,readout,st)
println(AC.energy_forces_virial(sys,calc;parallel=true).energy)

# Synthetic labels; training_rows accepts plain eV/Å numbers or Unitful labels.
# For ExtXYZ data, pass its energy, per-atom forces and virial properties here.
target_ps=A.install_readout(readout,1.01A.linear_parameters(readout))
target=A.ace_evaluate(input.positions,model,target_ps,input.environment,input.cell)
design=A.design_matrix(input.positions,model,readout,input.environment,input.cell)
labels=(energy=target.energy*u"eV",forces=target.forces*u"eV/Å",virial=target.virial*u"eV")
rows=A.training_rows(design,labels;weights=(energy=1.,forces=1.,virial=0.1))
coeff=A.fit_linear(rows.matrix,rows.labels;row_weights=rows.row_weights,
                   prior=A.smoothness_prior(model),regularization=1e-10)
fitted=A.install_readout(readout,coeff)
AC.set_parameters!(calc,fitted) # inference precontraction refreshes on next call
println("Linear prediction residual: ",maximum(abs,rows.matrix*coeff-rows.labels))

# Nonlinear training stays in original parameter coordinates. The Lux state is
# fixed and contains no reusable scratch. Topology is held fixed for this step.
layer=lux_layer(model;parameters=readout)
ps,lux_state=LuxCore.setup(LuxCore.Random.Xoshiro(42),layer)
loss(parameters)=sum(abs2,first(LuxCore.apply(layer,input,parameters,lux_state)).forces-target.forces)
gradient=Zygote.gradient(loss,ps)[1]
step(x::NamedTuple,g::NamedTuple,rate)=map((a,b) -> step(a,b,rate),x,g)
step(x::AbstractArray,g::AbstractArray,rate)=x-rate*g
updated=step(ps,gradient,1e-8)
println("Force loss before/after one nonlinear gradient step: ",(loss(ps),loss(updated)))
