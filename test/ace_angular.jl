@testset "explicit angular SIMD blocks and tiled complete forces" begin
    for L in (0,1,2,4,6), W in (1,2,4,8,16), T in (Float32,Float64)
        Rs=[SVector{3,T}(sin(i+0.3),cos(i*0.7),sin(2i+0.9)) for i=1:19]
        Y=zeros(T,length(Rs),(L+1)^2); dY=zeros(SVector{3,T},size(Y))
        block=ACECore.BlockAngular(L;width=W,precision=T)
        ACECore.angular_block!(Y,dY,block,Rs)
        ref=SpheriCart.SphericalHarmonics(max(1,L);normalisation=:L2,static=false,T)
        if L==0
            # Both SpheriCart 0.2.4 gradient kernels assume lmax>=1. Extract
            # the constant harmonic from lmax=1 and retain its zero derivative.
            full,dfull=SpheriCart.compute_with_gradients(ref,Rs)
            yr=full[:,1:1]; dyr=dfull[:,1:1]
        else
            yr,dyr=SpheriCart.compute_with_gradients(ref,Rs)
        end
        tol=T===Float64 ? 2e-12 : 2e-5
        @test Y ≈ yr atol=tol rtol=tol
        @test dY ≈ dyr atol=tol rtol=tol
    end
    model,ps,st=ACECore.load_resolved(ACE_MODEL_PATH)
    H=[4.7 0.5 0.1;0 4.8 0.3;0 0 4.6]
    sys=TestSystem([0.1 2.3 1.1;0.2 0.7 2.4;0.3 1.2 1.5],H,(true,true,true),[:Ti,:Al,:Ti])
    for width in (2,4,8), backend in (:sphericart,:generated,:block), fused in (false,true)
        inp=ACECore.ace_inputs(sys,model;environment=ACECore.FullEnvironment(;width))
        ref=ACECore.ace_evaluate(inp.positions,model,ps,inp.environment,inp.cell)
        w=ACECore.CPUExecutor(model,inp.environment;angular_backend=backend,slots=2)
        ACECore.execute_prepared!(w,inp.positions,inp.cell,ps;fused,parallel=true)
        y=ACECore.prepared_result(w)
        @test y.energy ≈ ref.energy atol=1e-10
        @test y.forces ≈ ref.forces atol=1e-10
        @test y.virial ≈ ref.virial atol=1e-10
    end
end
