using PotentialFrameWork, SIMD, InteractiveUtils, AtomsBase
const H=ChemicalSpecies(:H)
const P=bind_species(SmoothPair([H],ones(1,1),2.),[H])
const V=Vec{4,Float64}
const B=PairBlock(H,(H,H,H,H),V(1.),V(0.2),V(0.3),Vec{4,Bool}(true),Vec{4,Bool}(true))
# The callback is fully specialised: this output should contain vector fmul/fadd.
code_llvm(stdout,P,Tuple{typeof(B),EnergyForcesVirial};debuginfo=:none)
