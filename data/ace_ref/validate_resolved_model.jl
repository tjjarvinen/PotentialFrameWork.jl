# julia --project=data/ace_ref data/ace_ref/validate_resolved_model.jl
# Package under development is loaded alongside the separate pinned oracle.
using ACEpotentials # load pinned transitive dependencies before the core project
pushfirst!(LOAD_PATH,normpath(joinpath(@__DIR__,"../..")))
using PotentialFrameWork, ACEpotentials, ExtXYZ, AtomsCalculators, Unitful, LinearAlgebra, Test
const A=PotentialFrameWork.ACE
model,ps,st=A.load_resolved(joinpath(@__DIR__,"resolved_model_v1.toml"))
oracle,_=ACEpotentials.load_model(joinpath(@__DIR__,"reference_model.json"))
frames=ExtXYZ.load(joinpath(@__DIR__,"reference_results.xyz"))
errors=zeros(3)
@testset "resolved radial values and derivatives" begin
    for (family,reference) in ((model.radial,oracle.model.rbasis),(model.pair_radial,oracle.model.pairbasis)),
        i=1:2,j=1:2,r in range(0.1,5.49;length=71)
        a,da=ACEpotentials.Models.evaluate_ed(reference,r,model.atomic_numbers[i],model.atomic_numbers[j],NamedTuple(),NamedTuple())
        b,db=A.radial_ed(family[i,j],r)
        @test b ≈ a atol=1e-11 rtol=1e-11
        @test db ≈ da atol=1e-10 rtol=1e-10
    end
end
@testset "all 33 independent ACE predictions" begin
    for sys in frames
        input=A.ace_inputs(sys,model)
        r=A.ace_evaluate(input.positions,model,ps,input.environment,input.cell)
        ref=AtomsCalculators.energy_forces_virial(sys,oracle)
        e=ustrip(u"eV",ref.energy)
        f=hcat([ustrip.(u"eV/Å",v) for v in ref.forces]...)
        w=ustrip.(u"eV",ref.virial)
        errors .= max.(errors,[abs(r.energy-e),maximum(abs,r.forces-f),maximum(abs,r.virial-w)])
        @test r.energy ≈ e atol=1e-8 rtol=0
        @test r.forces ≈ f atol=1e-9 rtol=1e-10
        @test r.virial ≈ w atol=1e-8 rtol=1e-10
        @test sum(r.forces;dims=2) ≈ zeros(3,1) atol=1e-10
    end
end
println("Maximum independent E/F/W discrepancies: ",errors)

@testset "original-coordinate design and ACEfit solve" begin
    matrices=Matrix{Float64}[]; labels=Vector{Float64}[]
    for sys in frames[1:6]
        input=A.ace_inputs(sys,model)
        design=A.design_matrix(input.positions,model,ps,input.environment,input.cell)
        r=A.ace_evaluate(input.positions,model,ps,input.environment,input.cell)
        y=vcat(r.energy,vec(r.forces),vec(r.virial))
        @test design.matrix*A.linear_parameters(ps)+design.baseline ≈ y atol=1e-9 rtol=1e-10
        # Select energy, independent endpoint forces and diagonal virials. This
        # intentionally omits labels and leaves a rank-deficient linear problem.
        selected=[1,2,3,4,8,12,16]
        push!(matrices,design.matrix[selected,:]); push!(labels,(y-design.baseline)[selected])
    end
    matrix=vcat(matrices...); y=vcat(labels...)
    weights=1 .+ 0.2sin.(eachindex(y))
    coeff=A.fit_linear(matrix,y;row_weights=weights,regularization=1e-9,prior=I)
    @test matrix*coeff ≈ y atol=1e-6 rtol=1e-6
    toy=[1. 0 1;0 2 1;2 1 0;1 1 1]; labels=[1.,-2.,0.4,0.7]; prior=Diagonal([1.,2.,3.]); lambda=0.2
    fitted=A.fit_linear(toy,labels;prior,regularization=lambda)
    @test fitted ≈ (transpose(toy)*toy+lambda^2*transpose(prior)*prior)\(transpose(toy)*labels) atol=1e-12
    println("ACEfit ",pkgversion(ACEpotentials.ACEfit)," weighted E/F/W prediction recovery passed")
end
