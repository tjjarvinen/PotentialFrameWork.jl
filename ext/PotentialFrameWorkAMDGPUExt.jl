module PotentialFrameWorkAMDGPUExt
using PotentialFrameWork, AMDGPU
import PotentialFrameWork.ACE: storage_device, check_execution_device, automatic_gpu_strategy
storage_device(a::AMDGPU.ROCArray)=AMDGPU.device(a)
function check_execution_device(a::AMDGPU.ROCArray)
    AMDGPU.device(a)==AMDGPU.device() || throw(ArgumentError("select the adapted model's AMD device before preparing or executing GPU ACE"))
    nothing
end
# Calibrated only on the RX 6800 XT: TiAl at N=128/1024, both spline families.
# N=11664 showed no robust compact gain. Larger/custom models stay manual.
function automatic_gpu_strategy(m,reference::AMDGPU.ROCArray,n)
    small_splines=m.precision===Float64 && 128<=n<=1024 && m.source.lmax<=2 &&
        length(m.source.density)<=27 && length(m.source.products)<=167 &&
        all(ch->ch isa PotentialFrameWork.ACE.RadialChannel,m.source.radial) &&
        all(ch->ch isa PotentialFrameWork.ACE.RadialChannel,m.source.pair_radial)
    small_splines && occursin("AMD Radeon RX 6800 XT",string(AMDGPU.device(reference))) ? :compact : :staged
end
end
