const U=A.Unitful
const AB=A.AtomsBase
function audit_structure(z)
    data=bcc_geometry(4,z); X=data.positions
    H=ntuple(a->A.SVector{3}(data.cell[:,a])*U.u"Å",3)
    AB.atomic_system([AB.ChemicalSpecies(data.species[i])=>A.SVector{3}(X[:,i])*U.u"Å" for i=1:size(X,2)],H,(true,true,true))
end
