struct InvalidPrecisionUnits <: AtomsBase.AbstractSystem{3} end
Base.length(::InvalidPrecisionUnits)=1
AtomsBase.position(::InvalidPrecisionUnits,i::Integer)=fill(1.0u"fs",3)
AtomsBase.cell_vectors(::InvalidPrecisionUnits)=ntuple(j->SVector{3}(ntuple(a->a==j ? 16. : 0.,3))*u"Å",3)
AtomsBase.species(::InvalidPrecisionUnits,i::Integer)=ChemicalSpecies(:Ti)
AtomsBase.periodicity(::InvalidPrecisionUnits)=(false,false,false)
using SpheriCart
@testset "ACE precision adaptation and device primitives" begin
    A=PotentialFrameWork.ACE
    m,p,_=A.load_resolved(joinpath(@__DIR__,"../data/ace_ref/resolved_model_v1.toml"))
    b=A.adapt_model(m,p;precision=Float32)
    @test eltype(b.model.reference)===Float32
    @test eltype(b.model.coupling)===Float32
    @test eltype(b.parameters.WB)===Float32
    @test all(ch->eltype(ch.nodes)===Float32 && eltype(ch.coefficients)===Float32 && typeof(ch.transform.a)===Float32,b.model.radial)
    b.parameters.WB[1]=0; b.model.radial[1].nodes[1]=-.9
    @test p.WB[1]!=0 && m.radial[1].nodes[1]==-1
    back=A.adapt_model(A.adapt_model(m,p;precision=Float32).model,A.adapt_model(m,p;precision=Float32).parameters;precision=Float64)
    mixed=A.ACECalculator(A.adapt_model(m,p;precision=Float32).model,p;precision=Float32,angular_backend=:block)
    @test eltype(mixed.parameters.WB)===Float32
    @test eltype(mixed.executor.readout)===Float32
    @test back.model.coefficient_precision===Float32
    @test back.parameters.WB==Float64.(Float32.(p.WB))
    for T in (Float64,Float32), L in (0,1,2,4,6)
        Rs=[SVector{3,T}(r) for r in ((1.,0.,0.),(0.,0.,2.),(.7,-.9,1.3),(-1.2,.8,.04))]
        F=A.BlockAngular(L;precision=T,width=1).normalization
        ref,dref=A.angular_ed(A.angular_basis(L,T),Rs)
        for i in eachindex(Rs),k=1:(L+1)^2
            jet=A.gpu_angular(Rs[i],k,F)
            @test jet.v ≈ ref[i,k] atol=T===Float64 ? 2e-13 : 4e-6
            @test collect(jet.g) ≈ collect(dref[i,k]) atol=T===Float64 ? 2e-12 : 1e-5
        end
    end
    tables=A.pack_model(m,Array)
    for family in (:radial,:pair), iz=1:2,jz=1:2,r in (0.8,1.4,2.9,5.499999,5.5,5.6)
        ch=family===:radial ? m.radial[iz,jz] : m.pair_radial[iz,jz]
        values,derivatives=A.radial_ed(ch,r)
        for n in eachindex(values)
            v,d=A.gpu_radial(getproperty(tables,family),r,n,iz,jz)
            @test v ≈ values[n] atol=2e-12
            @test d ≈ derivatives[n] atol=2e-11
        end
    end
    # Invert the imported monotone transform to sample exact spline interval
    # boundaries and neighboring floating-point distances independently.
    ch=m.radial[1,1]
    for index in (2,div(length(ch.nodes),2),length(ch.nodes)-1)
        target=ch.nodes[index]; lo=ch.transform.rin; hi=ch.transform.rcut
        for _=1:70
            mid=(lo+hi)/2
            if A.transformed(ch.transform,mid)<target; lo=mid; else; hi=mid; end
        end
        for r in (prevfloat(hi),hi,nextfloat(hi))
            v,d=A.radial_ed(ch,r)
            for n in eachindex(v)
                value,derivative=A.gpu_radial(tables.radial,r,n,1,1)
                @test value ≈ v[n] atol=2e-12
                @test derivative ≈ d[n] atol=2e-11
            end
        end
    end
    sys=TestSystem([0. 2.8;0. .3;0. .2],Matrix{Float64}(I,3,3)*16,(true,true,true),[:Ti,:Al])
    report=A.compare_precision(m,p,sys)
    @test occursin("ACE precision comparison",sprint(show,MIME"text/plain"(),report))
    mktemp() do path,io
        close(io); cb=A.adapt_model(m,p;precision=Float32)
        calc=A.ACECalculator(cb.model,cb.parameters;precision=Float32,angular_backend=:block)
        @test eltype(calc.executor.readout)===Float32
        A.save_calculator(path,calc)
        loaded=A.ACECalculator(path;precision=Float64,angular_backend=:block)
        @test loaded.model.coefficient_precision===Float32
        @test A.compare_precision(loaded.model,loaded.parameters,sys).metadata.source_precision=="Float32"
    end
    @test length(report.structures)==1
    @test report.aggregate[:status]===:descriptive
    @test report.aggregate[:forces].components==6
    @test isfinite(report.aggregate[:forces].pooled_rms)
    @test A.compare_precision(m,p,[]).aggregate[:status]===:empty
    @test A.compare_precision(m,p,[sys,sys];tolerances=(energy=0.,forces=0.,virial=0.)).aggregate[:status]===:failed
    prep_failed=A.compare_precision(m,p,sys;array_type=x->error("unavailable test storage"))
    @test prep_failed.aggregate[:status]===:failed
    @test all(occursin("preparation failed",e) for e in prep_failed.structures[1].errors)
    invalid_units=A.compare_precision(m,p,[sys,InvalidPrecisionUnits(),sys])
    @test length(invalid_units.structures)==3
    @test invalid_units.structures[2].status===:failed
    @test isinf(invalid_units.aggregate[:forces].pooled_rms)
    bad=deepcopy(p); bad.WB[1]=NaN
    failed=A.compare_precision(m,bad,[sys,sys])
    @test failed.aggregate[:status]===:failed
    @test isinf(failed.aggregate[:forces].pooled_rms)
end
