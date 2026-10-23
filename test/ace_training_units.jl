@testset "ACE training label units" begin
    # Distinct components expose atom/component ordering and baseline mistakes.
    energy=7.25
    forces=reshape(collect(1.:6.),3,2)
    virial=reshape(collect(-4.:4.),3,3)
    design=(matrix=reshape(collect(1.:48.),16,3),baseline=collect(0.:15.))
    labels=(;energy,forces,virial)
    weights=(energy=2.,forces=0.3,virial=0.7)
    expected=ACECore.training_rows(design,labels;weights)
    @test expected.labels==vcat(energy,vec(forces),vec(virial))-design.baseline
    @test expected.row_weights==vcat(2.,fill(0.3,6),fill(0.7,9))
    for (eu,fu,vu) in ((u"eV",u"eV/Å",u"eV"),(u"J",u"N",u"J"),
                       (u"meV",u"eV/nm",u"meV"))
        unitful=(energy=uconvert(eu,energy*u"eV"),
                 forces=uconvert.(fu,forces*u"eV/Å"),
                 virial=uconvert.(vu,virial*u"eV"))
        saved=deepcopy(unitful)
        for f in (unitful.forces,[copy(v) for v in eachcol(unitful.forces)])
            rows=ACECore.training_rows(design,merge(unitful,(forces=f,));weights)
            @test rows.labels ≈ expected.labels
            @test rows.matrix==expected.matrix
            @test rows.row_weights==expected.row_weights
        end
        @test unitful==saved
    end
    mixed=(energy=energy*u"eV",forces=[forces[:,1]*u"eV/Å",forces[:,2]],virial=virial)
    @test ACECore.training_rows(design,mixed;weights).labels==expected.labels
    for key in (:energy,:forces,:virial)
        invalid=getproperty(labels,key)*u"s"
        @test_throws Unitful.DimensionError ACECore.training_rows(design,merge(labels,NamedTuple{(key,)}((invalid,))))
    end
    @test_throws DimensionMismatch ACECore.training_rows(design,(forces=[[1.,2.],[3.,4.]],))
    @test_throws DimensionMismatch ACECore.training_rows(design,(forces=zeros(3,3)*u"eV/Å",))
    selected=ACECore.training_rows(design,(energy=energy*u"eV",forces=nothing))
    @test selected.labels==[energy-design.baseline[1]]
    @test size(selected.matrix)==(1,3)
    @test isempty(ACECore.training_rows(design,(energy=nothing,)).labels)
    for fit_pair in (true,false)
        # AtomsBase geometry in nm; calculator outputs have Unitful quantities
        # and a vector of per-atom forces, matching the ExtXYZ boundary shape.
        model,ps,st=ACECore.load_resolved(ACE_MODEL_PATH)
        sys=AtomsBase.atomic_system([:Ti=>[0.,0.,0.]*u"nm",:Al=>[0.23,0.07,0.12]*u"nm"],
            ([1.,0.,0.]*u"nm",[0.,1.,0.]*u"nm",[0.,0.,1.]*u"nm"),(false,false,false))
        input=ACECore.ace_inputs(sys,model)
        @test input.positions[:,2] ≈ [2.3,0.7,1.2]
        d=ACECore.design_matrix(input.positions,model,ps,input.environment,input.cell;fit_pair)
        target=AC.energy_forces_virial(sys,ACECore.ACECalculator(model,ps,st))
        rows=ACECore.training_rows(d,target)
        coeff=fit_pair ? ACECore.linear_parameters(ps) : vec(ps.WB)
        @test rows.matrix*coeff ≈ rows.labels atol=1e-10
    end
    empty_design=(matrix=zeros(Float32,10,2),baseline=zeros(Float32,10))
    rows=ACECore.training_rows(empty_design,(energy=1f0*u"eV",forces=Vector{Float32}[]))
    @test rows.labels==Float32[1]
    @test eltype(rows.labels)==Float32
end
