module PotentialFrameWorkACEfitExt
using PotentialFrameWork, ACEfit, LinearAlgebra
import PotentialFrameWork.ACE: fit_linear

"""ACEfit matrix solve. `row_weights` multiply residuals directly (not squared
weights); omit rows for missing labels. Subtract fixed baseline before this call.
`regularization` multiplies `prior` in the augmented QR objective.
"""
function fit_linear(A::AbstractMatrix,y::AbstractVector;
                    row_weights=ones(length(y)),prior=I,regularization=0.,solver=nothing)
    size(A,1)==length(y)==length(row_weights) || throw(DimensionMismatch("fit rows"))
    all(w -> isfinite(w) && w>=0,row_weights) || throw(ArgumentError("invalid row weights"))
    s=isnothing(solver) ? ACEfit.QR(;lambda=regularization,P=prior) : solver
    weighted_matrix=A.*row_weights; weighted_labels=y.*row_weights
    result=if s isa ACEfit.QR && !iszero(s.lambda) && s.P isa AbstractMatrix && size(s.P,1)!=size(A,2)
        # ACEfit 0.3.1's QR appends ncoeff zero targets, assuming a square P.
        # R_old*T can be rectangular. Preserve its exact operator action and
        # use the actual ACEfit QR solve on the explicitly augmented objective.
        size(s.P,2)==size(A,2) || throw(DimensionMismatch("prior coefficient columns"))
        augmented=vcat(weighted_matrix,s.lambda*s.P)
        targets=vcat(weighted_labels,zeros(eltype(weighted_labels),size(s.P,1)))
        ACEfit.solve(ACEfit.QR(;lambda=0.,P=I),augmented,targets)
    else
        ACEfit.solve(s,weighted_matrix,weighted_labels)
    end
    result["C"]
end
end
