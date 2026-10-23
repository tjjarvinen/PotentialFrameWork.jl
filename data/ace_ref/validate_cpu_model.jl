using ACEpotentials
pushfirst!(LOAD_PATH,normpath(joinpath(@__DIR__,"../..")))
using PotentialFrameWork, ACEpotentials, ExtXYZ, Test
const A=PotentialFrameWork.ACE
model,ps,st=A.load_resolved(joinpath(@__DIR__,"resolved_model_v1.toml"))
frames=ExtXYZ.load(joinpath(@__DIR__,"reference_results.xyz"))
@testset "all-frame scalar/SIMD/angular/thread agreement" begin
    for width in (1,2,4,8,16),backend in (:sphericart,:generated,:block)
        env=A.FullEnvironment(;width)
        w=A.CPUExecutor(model,env;slots=4,angular_backend=backend)
        for sys in frames
            input=A.ace_inputs(sys,model;environment=env)
            reference=A.ace_evaluate(input.positions,model,ps,env,input.cell)
            A.prepare_executor!(w,env)
            for fused in (false,true),parallel in (false,true)
                A.execute_prepared!(w,input.positions,input.cell,ps;parallel,fused)
                result=A.prepared_result(w)
                @test result.energy ≈ reference.energy atol=1e-8 rtol=0
                @test result.forces ≈ reference.forces atol=1e-9 rtol=1e-10
                @test result.virial ≈ reference.virial atol=1e-8 rtol=1e-10
            end
        end
    end
end
println("Julia ",VERSION,", threads=",Threads.nthreads())
