using UnicodePlots
@testset "optional ACE terminal plots" begin
    A=PotentialFrameWork.ACE
    m,ps,st=A.load_resolved(ACE_MODEL_PATH)
    @test Base.get_extension(PotentialFrameWork,:PotentialFrameWorkUnicodePlotsExt)!==nothing
    for kind in (:orders,:weights,:coefficients,:coupling,:counts,:radial)
        @test A.plot_basis(m,ps;kind) isa UnicodePlots.Plot
    end
    empty=A.prune_basis(m,ps,A.BasisSelection(m;indices=Int[]))
    for kind in (:orders,:weights,:coefficients,:coupling,:counts,:radial)
        @test A.plot_basis(empty.model,empty.parameters;kind) isa UnicodePlots.Plot
    end
    samples=A.radial_samples(m;channels=[1,3],points=11)
    @test size(samples.values)==(11,2)
    @test size(A.radial_samples(empty.model;family=:pair_radial,points=11).values)==(11,6)
    @test samples.values[4,:] ≈ A.radial_values(m.radial[1,1],samples.distances[4])[[1,3]]
    @test_throws ArgumentError A.plot_basis(m;kind=:coefficients)
end
