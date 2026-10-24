using ACEfit, LinearAlgebra, Unitful, AtomsBase
include("ace_interface_plotting.jl")
@testset "selected ACEfit regression" begin
    A=PotentialFrameWork.ACE
    m,ps,st=A.load_resolved(ACE_MODEL_PATH)
    reduced=A.prune_basis(m,ps,A.BasisSelection(m;indices=[1,4,9,17]))
    model=A.project_species(reduced.model,[:Al],0:2)
    parameters=reduced.parameters
    coefficients=A.linear_parameters(parameters).*1.07
    target_parameters=A.install_readout(parameters,coefficients)
    matrices=Matrix{Float64}[]; targets=Vector{Float64}[]; weights=Vector{Float64}[]
    for distance in (1.9,2.2,2.7,3.1)
        sys=AtomsBase.isolated_system([:Ti=>[0.,0.,0.]*u"Å",:Al=>[distance,0.3,0.7]*u"Å",:Ti=>[0.7,2.1,0.4]*u"Å"])
        input=A.ace_inputs(sys,model)
        d=A.design_matrix(input.positions,model,parameters,input.environment,input.cell;centers=[1,3])
        target=A.ace_evaluate(input.positions,model,target_parameters,input.environment,input.cell;centers=[1,3])
        rows=A.training_rows(d,(energy=target.energy*u"eV",forces=target.forces*u"eV/Å",virial=target.virial*u"eV");
            weights=(energy=2.,forces=0.7,virial=0.3))
        push!(matrices,rows.matrix); push!(targets,rows.labels); push!(weights,rows.row_weights)
    end
    design=vcat(matrices...); labels=vcat(targets...); row_weights=vcat(weights...)
    fitted=A.fit_linear(design,labels;row_weights,regularization=1e-10,prior=I)
    @test design*fitted ≈ labels atol=1e-7 rtol=1e-7
    @test A.linear_parameters(A.install_readout(parameters,fitted))==fitted
    # A genuine non-diagonal coordinate map makes R_old*T rectangular. Check
    # the actual solver objective, not just an algebraic map identity.
    original_prior=[1. 0.2 0;0.1 2 0.3;0 0.4 3]
    transform=[1. 0.2;0.3 1;0.1 -0.4]
    prior=A.mapped_prior(original_prior,transform)
    toy=[1. 0.2;0.3 2;1.2 -0.5;0.7 0.9]; target=[1.,-2.,0.4,0.7]; lambda=0.2
    expected=(transpose(toy)*toy+lambda^2*transpose(prior)*prior)\(transpose(toy)*target)
    @test A.fit_linear(toy,target;prior,regularization=lambda) ≈ expected atol=1e-12
    @test A.fit_linear(toy,target;solver=ACEfit.QR(;lambda,P=prior)) ≈ expected atol=1e-12
end
