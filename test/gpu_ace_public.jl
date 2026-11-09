# Explicit external data: julia --project=<app> test/gpu_ace_public.jl DATA_DIR [--amdgpu]
using PotentialFrameWork, Test, TOML, LinearAlgebra
const PA=PotentialFrameWork.ACE
const PN=PotentialFrameWork.NeighbourList
const AB=PA.AtomsBase
const U=PA.Unitful
storage=if "--amdgpu" in ARGS
    @eval using AMDGPU
    AMDGPU.functional() || error("requested AMD unavailable")
    AMDGPU.device!(only(filter(d->occursin("6800 XT",string(d)),AMDGPU.devices())))
    AMDGPU.ROCArray
else
    Array
end
reports=Dict[]
@testset "independent public ACE models" begin
for name in ("ti","tialv")
    dir=ARGS[1]; model,ps,_=PA.load_resolved(joinpath(dir,name*"-resolved.toml"))
    refs=TOML.parsefile(joinpath(dir,name*"-oracle.toml"))
    @test refs["source_sha256"]==model.source_sha256
    @test model.coordinate_kind===:legacy_PI
    @test_throws ArgumentError PA.descriptors(model)
    for (frame,record) in enumerate(refs["records"])
        n=length(record["species"]); X=reshape(record["positions"],3,n); H=reshape(record["cell"],3,3)
        cells=ntuple(j->PA.SVector{3}(H[:,j])*U.u"Å",3)
        sys=AB.atomic_system([AB.ChemicalSpecies(record["species"][i])=>PA.SVector{3}(X[:,i])*U.u"Å" for i=1:n],cells,(true,true,true))
        ref=(;energy=record["energy"],forces=reshape(record["forces"],3,n),virial=reshape(record["virial"],3,3))
        for T in (Float64,Float32)
            bundle=PA.adapt_model(model,ps;precision=T,array_type=storage)
            result=if storage===Array
                PA.precision_evaluate(sys,bundle,(PA.AC.Energy(),PA.AC.Forces(),PA.AC.Virial()))
            else
                PA.gpu_result(PA.execute_gpu!(PA.GPUExecutor(sys,bundle.model,bundle.parameters;batch_size=16)))
            end
            e=result.energy isa Number ? result.energy : only(result.energy)
            ef=maximum(abs,result.forces-ref.forces); ew=maximum(abs,result.virial-ref.virial); ee=abs(e-ref.energy)
            println(name," frame=",frame," ",T," E/F/W=",(ee,ef,ew))
            cpu=PA.adapt_model(model,ps;precision=T)
            rounded=PA.precision_evaluate(sys,cpu,(PA.AC.Energy(),PA.AC.Forces(),PA.AC.Virial()))
            implementation_force=maximum(abs,result.forces-rounded.forces;init=0.)
            implementation_virial=maximum(abs,result.virial-rounded.virial;init=0.)
            if T===Float64
                @test implementation_force<=2e-9
                @test implementation_virial<=1e-8
                @test e ≈ ref.energy atol=2e-7 rtol=1e-10
                @test result.forces ≈ ref.forces atol=2e-7 rtol=1e-9
                @test result.virial ≈ ref.virial atol=2e-6 rtol=1e-9
            else
                # Application precision sensitivity is descriptive here; compare
                # the device separately to CPU on the same rounded model/input.
                @test result.forces ≈ rounded.forces atol=2e-3 rtol=1e-3
                @test result.virial ≈ rounded.virial atol=2e-2 rtol=1e-3
            end
            push!(reports,Dict("model"=>name,"frame"=>frame,"atoms"=>n,"precision"=>string(T),"energy_abs_eV"=>ee,"force_max_eV_per_A"=>ef,"virial_max_eV"=>ew,"implementation_force_max_eV_per_A"=>implementation_force,"implementation_virial_max_eV"=>implementation_virial))
        end
    end
end
end
if haskey(ENV,"PFW_PUBLIC_REPORT")
    open(ENV["PFW_PUBLIC_REPORT"],"w") do io; TOML.print(io,Dict("cases"=>reports,"storage"=>string(storage),"device"=>storage===Array ? "CPU" : string(AMDGPU.device()),"julia"=>string(VERSION))); end
end

if storage===Array
    hashes=Dict(name=>bytes2hex(PA.SHA.sha256(read(joinpath(ARGS[1],name*"-resolved.toml")))) for name in ("ti","tialv"))
    open(joinpath(ARGS[1],"verification-cpu.toml"),"w") do io
        TOML.print(io,Dict("schema"=>"PotentialFrameWork.PublicACEVerification","version"=>1,"resolved_sha256"=>hashes,"julia"=>string(VERSION)))
    end
end
