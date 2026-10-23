using TOML, SpheriCart, ForwardDiff
import AtomsCalculators as AC

@testset "ACE calculator file interface and trained checkpoints" begin
    A=PotentialFrameWork.ACE
    resolved=joinpath(@__DIR__,"../data/ace_ref/resolved_model_v1.toml")
    recipe=joinpath(dirname(resolved),"reference_model.json")
    original=A.ACECalculator(resolved;slots=1)
    from_recipe=A.ACECalculator(recipe;resolved_path=resolved,slots=1)
    @test AC.get_parameters(from_recipe)==AC.get_parameters(original)
    H=[4.7 0.5 0.1;0 4.8 0.3;0 0 4.6]
    sys=TestSystem([0.1 2.3 1.1;0.2 0.7 2.4;0.3 1.2 1.5],H,(true,true,true),[:Ti,:Al,:Ti])
    before=AC.energy_forces_virial(sys,original) # populate caches before tuning
    trained=A.trainable_parameters(original.model,AC.get_parameters(original))
    trained.WB .*= 0.83; trained.Wpair .*= 1.12
    trained.radial.scale[:,1,2] .*= 0.91
    trained.pair.scale[1,2,1]=0.0
    trained.pair.scale[2,1,2]=1.17
    state=(epoch=Int32(7),rate=Float32(0.002),tag="tuned",extra=(nothing,true,zeros(Float32,0,2)))
    expected=AC.calculate((AC.Energy(),AC.Forces(),AC.Virial()),sys,original;parameters=trained,state)
    @test expected.energy!=before.energy
    mktempdir() do dir
        path=joinpath(dir,"trained.toml")
        @test A.save_calculator(path,original;parameters=trained,state)==path
        @test original.parameters!=trained # explicit override does not install parameters
        @test original.state==NamedTuple()
        saved=TOML.parsefile(path)
        @test Set(keys(saved))==Set(["schema","schema_version","writer_julia","model","parameters","state"])
        @test !haskey(saved["model"],"parameters") # one authoritative tree
        # Load elsewhere, with no recipe/export sibling, and different scratch layout.
        relocated=joinpath(dir,"moved","checkpoint.ace")
        mkpath(dirname(relocated)); mv(path,relocated)
        restored=A.ACECalculator(relocated;width=8,slots=1,angular_backend=:block)
        @test restored.parameters==trained
        @test restored.state==state
        @test typeof(restored.state)==typeof(state)
        @test restored.parameters.WB!==trained.WB
        @test restored.model.source_sha256==original.model.source_sha256
        @test restored.model.coupling==original.model.coupling
        @test restored.model.products==original.model.products
        @test isempty(restored.environment.neighbours)
        actual=AC.energy_forces_virial(sys,restored)
        @test actual.energy ≈ expected.energy atol=1e-10u"eV"
        @test ustrip.(hcat(actual.forces...)) ≈ ustrip.(hcat(expected.forces...)) atol=1e-10
        @test actual.virial ≈ expected.virial atol=1e-10u"eV"
        # A loaded checkpoint remains usable with explicit AtomsCalculators AD.
        force_loss(t,c,ps)=sum(abs2,ustrip.(hcat(AC.calculate(AC.Forces(),sys,c;
            parameters=merge(ps,(WB=ps.WB.*t,))).forces...)))
        @test ForwardDiff.derivative(t -> force_loss(t,restored,restored.parameters),1.0) ≈
              ForwardDiff.derivative(t -> force_loss(t,original,trained),1.0) rtol=1e-12
        # Defaults save installed values, even after old executor weights were used.
        AC.set_parameters!(restored,merge(trained,(WB=trained.WB.*0.7,)))
        A.save_calculator(path,restored)
        twice=A.ACECalculator(path;slots=1)
        @test twice.parameters==AC.get_parameters(restored)
        a=AC.energy_forces_virial(sys,twice); b=AC.energy_forces_virial(sys,restored)
        @test a.energy ≈ b.energy atol=1e-10u"eV"
        @test ustrip.(hcat(a.forces...)) ≈ ustrip.(hcat(b.forces...)) atol=1e-10
        @test a.virial ≈ b.virial atol=1e-10u"eV"
        # Float32 parameter blocks retain their element type and multidimensional shapes.
        ps32=(WB=Float32.(trained.WB),Wpair=Float32.(trained.Wpair),
              radial=(scale=Float32.(trained.radial.scale),),pair=(scale=Float32.(trained.pair.scale),))
        A.save_calculator(path,original;parameters=ps32)
        single=A.ACECalculator(path;slots=1,precision=Float32)
        @test typeof(single.parameters)==typeof(ps32)
        @test single.parameters==ps32
        # Failed saves preserve an existing checkpoint; no silent dropped fields/tracers.
        intact=read(path)
        for (ps,st) in ((merge(trained,(WB=fill(NaN,size(trained.WB)),)),NamedTuple()),
                        (trained,(bad=Ref(1),)))
            @test_throws ArgumentError A.save_calculator(path,original;parameters=ps,state=st)
            @test read(path)==intact
        end
        @test_throws DimensionMismatch A.save_calculator(path,original;parameters=merge(trained,(radial=(scale=ones(2,2),),)))
        @test read(path)==intact
        @test_throws DimensionMismatch A.save_calculator(path,original;parameters=(WB=zeros(1,1),Wpair=trained.Wpair))
        dual=ForwardDiff.Dual.(trained.WB,1.0)
        @test_throws ArgumentError A.save_calculator(path,original;parameters=merge(trained,(WB=dual,)))
        @test read(path)==intact
        # The file content identifies checkpoints even with a recipe-like extension.
        jsonpath=joinpath(dir,"checkpoint.json")
        A.save_calculator(jsonpath,original)
        @test A.ACECalculator(jsonpath;slots=1).parameters==original.parameters
        # Explicit schema/shape errors fail before numerical execution.
        saved["schema_version"]=2
        open(io -> TOML.print(io,saved),path,"w")
        @test_throws ArgumentError A.ACECalculator(path)
        saved["schema_version"]=1
        saved["parameters"]["items"][1]["shape"]=[1,1]
        open(io -> TOML.print(io,saved),path,"w")
        @test_throws DimensionMismatch A.ACECalculator(path)
    end
end
