# Run in an environment containing PotentialFrameWork, ForwardDiff and Zygote.
using PotentialFrameWork, AtomsBase, Unitful, ForwardDiff, Zygote
import AtomsCalculators as AC

sys=AtomsBase.FlexibleSystem([
    AtomsBase.Atom(:Ar,[0.,0.,0.]*u"Å"),
    AtomsBase.Atom(:Ar,[1.,0.,0.]*u"Å")],
    [[10.,0.,0.],[0.,10.,0.],[0.,0.,10.]].*u"Å",(false,false,false))
p=SmoothPair([:Ar],reshape([0.1],1,1),2.)
parameters=AC.get_parameters(p) # (amplitudes = [0.1],)
energy(theta)=ustrip(u"eV",AC.calculate(AC.Energy(),sys,p;
    parameters=(amplitudes=theta,)).energy)
println("energy (eV): ",energy(parameters.amplitudes))
println("ForwardDiff: ",ForwardDiff.gradient(energy,parameters.amplitudes))
println("ChainRules via Zygote: ",Zygote.gradient(energy,parameters.amplitudes)[1])
