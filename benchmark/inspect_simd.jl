using PotentialFrameWork, SIMD, InteractiveUtils
using PotentialFrameWork.PairPotentials: BoundSmooth
const P=BoundSmooth(ones(2,2),2.)
const V=Vec{4,Float64}
const B=PairBlock(1,(1,1,1,1),V(1.),V(0.2),V(0.3),Vec{4,Bool}(true),Vec{4,Bool}(true))
# The callback is fully specialised: this output should contain vector fmul/fadd.
code_llvm(stdout,P,Tuple{typeof(B),EnergyForcesVirial};debuginfo=:none)
