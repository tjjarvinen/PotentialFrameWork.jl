# Run in an application environment with AMDGPU installed; core CPU use needs no
# AMDGPU installation. See docs/gpu_pairs.md for environment setup.
using PotentialFrameWork, AMDGPU, AtomsBuilder, AtomsSystems, Unitful
import AtomsCalculators as AC
const NL=PotentialFrameWork.NeighbourList
const PP=PotentialFrameWork.PairPotentials
AMDGPU.functional() || error("AMDGPU is unavailable")
AMDGPU.device!(only(filter(d->occursin("6800 XT",string(d)),AMDGPU.devices())))
backend=NL.KA.get_backend(ROCArray(zeros(1)))
sys=AtomsSystems.generic_system(repeat(AtomsBuilder.bulk(:Ar;cubic=true),(6,6,6)))
p=PP.lennard_jones()
for strategy in (:central,:candidate)
    config=NL.GPUConfig(;strategy,workgroup=64,tile=4) # skin=0: fresh cells
    input=PP.prepare_gpu_execution(sys,p,backend;precision=Float32,config)
    result=PP.execute_gpu!(input;request=(AC.Energy(),AC.Forces(),AC.Virial()))
    NL.KA.synchronize(backend)
    println(strategy," device energy = ",Array(result.energy))
    # Explicit parameters are refreshed each call; positions remain on device.
    updated=PP.execute_gpu!(input;parameters=(;epsilon=.012,sigma=3.4))
    NL.KA.synchronize(backend)
    println("updated device energy = ",Array(updated.energy))
    owning=PP.evaluate_gpu!(input.workspace,sys,p;request=(AC.Forces(),))
    println("unitful original-order force size = ",size(owning.forces))
end
# Buffered connectivity is opt-in. Updating input.positions keeps original order.
input=PP.prepare_gpu_execution(sys,p,backend;precision=Float32,config=NL.GPUConfig(;skin=.5))
PP.execute_gpu!(input)
input.positions .+= 0.001f0
PP.execute_gpu!(input)
NL.KA.synchronize(backend)
println("skin rebuilds = ",input.workspace.rebuilds)
