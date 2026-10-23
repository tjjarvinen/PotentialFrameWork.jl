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
    result=ACEfit.solve(s,A.*row_weights,y.*row_weights)
    result["C"]
end
end
