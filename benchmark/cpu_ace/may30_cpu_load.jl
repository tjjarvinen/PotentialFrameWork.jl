using PotentialFrameWork
const A=PotentialFrameWork.ACE
@assert isnothing(Base.get_extension(PotentialFrameWork,:PotentialFrameWorkAMDGPUExt))
@assert isnothing(Base.get_extension(PotentialFrameWork,:PotentialFrameWorkSpheriCartExt))
m,p,_=A.load_resolved(joinpath(pwd(),"data/may30/resolved.toml"))
d=A.TOML.parsefile(joinpath(pwd(),"data/may30/bcc-128.toml"))
X=reshape(d["positions"],3,:); H=reshape(d["cell"],3,3)
sys=A.AtomsBase.atomic_system([A.ChemicalSpecies(d["species"][i])=>A.SVector{3}(X[:,i])*A.u"Å" for i in axes(X,2)],ntuple(j->A.SVector{3}(H[:,j])*A.u"Å",3),(true,true,true))
c=A.ACECalculator(m,p;slots=1,angular_backend=:block,product_algorithm=:dag)
y=A.AC.energy_forces_virial(sys,c)
@assert abs(A.ustrip(A.u"eV",y.energy)-d["energy"])<2e-8
println("CPU-only loading and solid block E/F/V passed; no SpheriCart or vendor extension loaded")
