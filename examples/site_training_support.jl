# Small generated CPU systems shared by the two site-training examples.
module SiteTrainingSupport
using PotentialFrameWork, AtomsBase, Unitful, LinearAlgebra
const A=PotentialFrameWork.ACE
const NL=PotentialFrameWork.NeighbourList
struct TrainingSystem <: AtomsBase.AbstractSystem{3}
    X::Matrix{Float64}
    H::Matrix{Float64}
    periodic::NTuple{3,Bool}
    labels::Vector{ChemicalSpecies}
end
Base.length(s::TrainingSystem)=size(s.X,2)
AtomsBase.position(s::TrainingSystem,i::Integer)=s.X[:,i]*u"Å"
AtomsBase.cell_vectors(s::TrainingSystem)=ntuple(i -> s.H[:,i]*u"Å",3)
AtomsBase.periodicity(s::TrainingSystem)=s.periodic
AtomsBase.species(s::TrainingSystem,i::Integer)=s.labels[i]
function setup(;periodic=false)
    ch=A.CallableRadial(A.SmoothExponential(),A.RadialMetadata([(1,0),(1,1),(2,1),(1,2)];cutoff=3.2,ad=true);
        config=(cutoff=3.2,powers=[0,0,1,0]))
    l=A.ACEFeatureLayer(species=[:Si,:O],radial=ch,correlation=2,lmax=2,totaldegree=4,
        output_irreps=((L=0,parity=:even,channels=2),(L=1,parity=:odd,channels=2),
                       (L=1,parity=:even,channels=1),(L=2,parity=:even,channels=2)))
    ps=A.feature_parameters(l;custom_radial=(rate=[0.4],amplitude=[1.3]))
    for (b,W) in enumerate(values(ps.mixing)); W .= reshape(0.3sin.(b .+ (1:length(W))),size(W)); end
    sys=TrainingSystem([0.1 1.4 2.2;0.2 0.4 1.1;0.3 0.7 0.2],Matrix(4.0I,3,3),(periodic,periodic,periodic),ChemicalSpecies.([:Si,:O,:Si]))
    n=NL.prepare_neighbors(sys,3.2;channels=[1,2,1],species_count=2,width=4)
    l,ps,n
end
step(ps::NamedTuple,g,eta)=NamedTuple{keys(ps)}(map(k -> step(getproperty(ps,k),getproperty(g,k),eta),keys(ps)))
step(ps::AbstractArray,g,eta)=ps.-eta.*g
step(ps::Real,g,eta)=ps-eta*g
struct SmoothInvariantHead end
function (::SmoothInvariantHead)(h,ps)
    q=[sum(abs2,block) for block in h]
    u=ps.bias+sum(ps.weights.*q)
    ps.amplitude*tanh(u)+ps.quadratic*sum(q)^2
end

end
