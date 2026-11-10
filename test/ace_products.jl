@testset "prepared product algorithm: division-free values and weighted adjoints" begin
    A=ACECore
    products=[[1],[2,2],[1,2,3],[4,1,4,2],[2,2,2,3,4,1]]
    for T in (Float32,Float64)
        tol=T===Float32 ? 2e-5 : 3e-12
        for U0 in ([.7,-.3,1.2,-.6],[0.,-.3,1.2,.6],[0.,0.,1.2,-.6],[1e-3,-3.,.2,-.6]),
            V0 in (zeros(4),[.2,.3,-.4,.1]), pass in 1:2
            U=T.(U0).*T(pass); V=T.(V0)
            weights=T.([.2,0.,-.7,1.3,-.4])./T(pass)
            fwd=zeros(T,7,7); rev=similar(fwd)
            P=zeros(T,5); dU=zeros(T,4); dV=similar(dU)
            A.product_projection!(P,dU,dV,U,V,products,weights,nothing,fwd,rev)
            expected=(copy(P),copy(dU),copy(dV))
            # Dirty tables and adjoints must not survive repeated calls.
            fill!(fwd,42); fill!(rev,-17); fill!(dU,99); fill!(dV,-99)
            A.prepared_products!(P,dU,dV,U,V,products,weights,nothing,fwd,rev,:optimized)
            for (actual,ref) in zip((P,dU,dV),expected)
                @test actual ≈ ref atol=tol rtol=tol
            end
            # Independent direct expression, including repeated destinations.
            loss(u,v)=sum(weights[j]*prod(u[a]+v[a] for a in ids) for (j,ids) in enumerate(products))
            h=T===Float32 ? T(.002) : T(1e-5)
            # Central differences incur O(h²) truncation and O(eps/h) rounding.
            fd_tol=T===Float32 ? 2e-4 : 5e-8
            for a in (1,2,4)
                du=zeros(T,4); du[a]=h
                @test dU[a] ≈ (loss(U+du,V)-loss(U-du,V))/(2h) atol=fd_tol rtol=fd_tol
                @test dV[a] ≈ (loss(U,V+du)-loss(U,V-du))/(2h) atol=fd_tol rtol=fd_tol
            end
            for orders in (0:0,1:2,0:6)
                A.product_projection!(P,dU,dV,U,V,products,weights,orders,fwd,rev)
                expected=(copy(P),copy(dU),copy(dV))
                A.prepared_products!(P,dU,dV,U,V,products,weights,orders,fwd,rev,:optimized)
                @test (P,dU,dV)==expected
            end
        end
    end
    # Empty product sets occur after structural pruning; empty factors do not
    # occur in validated model products and are deliberately not introduced.
    for T in (Float32,Float64), algorithm in (:reference,:optimized,:auto)
        P=T[]; dU=ones(T,2); dV=ones(T,2)
        A.prepared_products!(P,dU,dV,T[1,2],T[3,4],Vector{Int}[],T[],nothing,
            zeros(T,1,1),zeros(T,1,1),algorithm)
        @test isempty(P) && all(iszero,dU) && all(iszero,dV)
    end
    m,ps,_=A.load_resolved(ACE_MODEL_PATH)
    sys=TestSystem([.1 2.3 1.1;.2 .7 2.4;.3 1.2 1.5],
        [3.7 .7 .1;0 3.8 .3;0 0 3.6],(true,true,true),[:Ti,:Al,:Ti])
    masked=A.mask_basis(m,A.BasisSelection(m;indices=1:2:size(m.coupling,1)))
    split=A.model_with(m;policy=A.policy_with(m.policy;species=(1,)))
    for model in (m,masked,split,A.project_species(m,[:Ti],1:2)), T in (Float32,Float64)
        b=A.adapt_model(model,ps;precision=T)
        input=A.ace_inputs(sys,b.model;environment=A.FullEnvironment(width=4))
        w=A.CPUExecutor(b.model,input.environment;precision=T,slots=2)
        for params in (b.parameters,A.install_readout(b.parameters,.7A.linear_parameters(b.parameters))), centers in (nothing,[2])
            A.execute_prepared!(w,input.positions,input.cell,params;product_algorithm=:reference,centers)
            ref=A.prepared_result(w)
            for algorithm in (:optimized,:auto)
                A.execute_prepared!(w,input.positions,input.cell,params;product_algorithm=algorithm,parallel=true,centers)
                r=A.prepared_result(w); tol=T===Float32 ? 3e-4 : 2e-10
                @test r.energy ≈ ref.energy atol=tol rtol=tol
                @test r.forces ≈ ref.forces atol=tol rtol=tol
                @test r.virial ≈ ref.virial atol=tol rtol=tol
            end
        end
        @test_throws ArgumentError A.execute_prepared!(w,input.positions,input.cell,b.parameters;product_algorithm=:wrong)
    end
    @test_throws ArgumentError A.CPUExecutor(m,A.FullEnvironment();product_algorithm=:wrong)
end

@testset "optimized prepared force and strain finite differences" begin
    A=ACECore
    m,ps,_=A.load_resolved(ACE_MODEL_PATH); m=A.components(m;reference=false)
    X=[.4 2.6 1.1;0. .2 2.8;0. .3 .1]
    H=[16. .4 .1;0 17. .3;0 0 18.]
    for periodic in (false,true)
        sys=TestSystem(X,H,(periodic,periodic,periodic),[:Ti,:Al,:Ti])
        input=A.ace_inputs(sys,m;environment=A.FullEnvironment(width=4))
        w=A.CPUExecutor(m,input.environment;product_algorithm=:optimized)
        evaluate(x,h)=begin
            A.execute_prepared!(w,x,h,ps)
            A.prepared_result(w)
        end
        y=evaluate(input.positions,input.cell); h=1e-5
        for a=1:3,i=1:3
            d=zeros(3,3); d[a,i]=h
            fd=-(evaluate(input.positions+d,input.cell).energy-evaluate(input.positions-d,input.cell).energy)/(2h)
            @test fd ≈ y.forces[a,i] atol=2e-7 rtol=2e-7
        end
        for a=1:3,b=1:3
            d=zeros(3,3); d[b,a]=h
            fd=-(evaluate((I+d)*input.positions,(I+d)*input.cell).energy-
                 evaluate((I-d)*input.positions,(I-d)*input.cell).energy)/(2h)
            @test fd ≈ y.virial[a,b] atol=2e-7 rtol=2e-7
        end
        calc=A.ACECalculator(m,ps;product_algorithm=:reference)
        r=A.AC.energy_forces_virial(sys,calc;product_algorithm=:optimized)
        @test ustrip(u"eV",r.energy) ≈ y.energy atol=2e-10
        @test reduce(hcat,[ustrip.(u"eV/Å",v) for v in r.forces]) ≈ y.forces atol=2e-10
        @test ustrip.(u"eV",r.virial) ≈ y.virial atol=2e-10
        p=A.NI.prepare_potential(m,[AtomsBase.species(sys,i) for i=1:length(sys)];width=4,product_algorithm=:reference)
        @test p.product_algorithm===:reference
    end
end
