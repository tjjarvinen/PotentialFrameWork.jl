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
            for algorithm in (:optimized,:auto,:dag)
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
    for periodic in (false,true), algorithm in (:optimized,:dag)
        sys=TestSystem(X,H,(periodic,periodic,periodic),[:Ti,:Al,:Ti])
        input=A.ace_inputs(sys,m;environment=A.FullEnvironment(width=4))
        w=A.CPUExecutor(m,input.environment;product_algorithm=algorithm)
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
        r=A.AC.energy_forces_virial(sys,calc;product_algorithm=algorithm)
        @test ustrip(u"eV",r.energy) ≈ y.energy atol=2e-10
        @test reduce(hcat,[ustrip.(u"eV/Å",v) for v in r.forces]) ≈ y.forces atol=2e-10
        @test ustrip.(u"eV",r.virial) ≈ y.virial atol=2e-10
        p=A.NI.prepare_potential(m,[AtomsBase.species(sys,i) for i=1:length(sys)];width=4,product_algorithm=:reference)
        @test p.product_algorithm===:reference
    end
end

@testset "product DAG values, weighted reverse and structural lifetime" begin
    A=ACECore
    # Duplicate terminals, canonical reordering, outputs used internally, equal
    # operands, repeated and unused leaves, every order through six.
    products=[[2],[2,1],[1,2],[2,2],[2,1,3],[1,2,3,4],
              [4,3,2,1,4],[2,1,2,3,4,4],[1,1],[1,1,1]]
    plan=A.ProductDAG(products,5)
    @test plan.terminals[2]==plan.terminals[3]
    @test A.ProductDAG(products,5).left==plan.left
    @test_throws ArgumentError A.ProductDAG([Int[]],5)
    @test_throws ArgumentError A.ProductDAG([[0]],5)
    for T in (Float32,Float64)
        tol=T===Float32 ? 2e-5 : 3e-12
        P=zeros(T,length(products)); dU=zeros(T,5); dV=similar(dU)
        fwd=zeros(T,7,7); rev=similar(fwd)
        values=zeros(T,5+length(plan.left)); adjoints=similar(values)
        for u in ([.7,-.3,1.2,-.6,3.],[0.,-.3,1.2,.6,0.],[0.,0.,1.2,-.6,2.]),
            v in (zeros(5),[.2,.3,-.4,.1,0.]), pass in 1:3
            U=T.(u).*T(pass); V=T.(v)
            weights=T[isodd(j+pass) ? 0 : sin(j)/pass for j in eachindex(products)]
            expected=Tuple[]
            for algorithm in (:reference,:optimized)
                A.prepared_products!(P,dU,dV,U,V,products,weights,nothing,fwd,rev,algorithm)
                push!(expected,(copy(P),copy(dU),copy(dV)))
            end
            fill!(values,T(NaN)); fill!(adjoints,T(NaN))
            A.prepared_products!(P,dU,dV,U,V,products,weights,nothing,fwd,rev,:dag,plan,values,adjoints)
            for reference in expected, (y,r) in zip((P,dU,dV),reference)
                @test y ≈ r atol=tol rtol=tol
            end
            @test P ≈ [prod(U[a]+V[a] for a in ids) for ids in products] atol=tol rtol=tol
            @test iszero(dU[5]) && iszero(dV[5])
            loss(u,v)=sum(weights[j]*prod(u[a]+v[a] for a in ids) for (j,ids) in enumerate(products))
            h=T===Float32 ? T(.002) : T(1e-5)
            ftol=T===Float32 ? 2e-4 : 5e-8
            for a in (1,2,4)
                delta=zeros(T,5); delta[a]=h
                @test dU[a] ≈ (loss(U+delta,V)-loss(U-delta,V))/(2h) atol=ftol rtol=ftol
                @test dV[a] ≈ (loss(U,V+delta)-loss(U,V-delta))/(2h) atol=ftol rtol=ftol
            end
            for orders in (0:0,1:2,0:6)
                A.product_projection!(P,dU,dV,U,V,products,weights,orders,fwd,rev)
                reference=(copy(P),copy(dU),copy(dV))
                A.prepared_products!(P,dU,dV,U,V,products,weights,orders,fwd,rev,:dag,plan,values,adjoints)
                @test (P,dU,dV)==reference
                @test A.effective_product_algorithm(:dag,orders)===:reference
            end
        end
        emptyplan=A.ProductDAG(Vector{Int}[],5)
        A.dag_products!(T[],dU,dV,ones(T,5),ones(T,5),T[],emptyplan,zeros(T,5),ones(T,5))
        @test all(iszero,dU) && all(iszero,dV)
    end
    m,ps,_=A.load_resolved(ACE_MODEL_PATH)
    sys=TestSystem([.1 2.3 1.1;.2 .7 2.4;.3 1.2 1.5],
        [3.7 .7 .1;0 3.8 .3;0 0 3.6],(true,true,true),[:Ti,:Al,:Ti])
    input=A.ace_inputs(sys,m;environment=A.FullEnvironment(width=4))
    w=A.CPUExecutor(m,input.environment;product_algorithm=:dag,slots=2)
    p=w.prepared; plan=p.product_plan
    @test all(s->length(s.dag_values)==plan.nleaves+length(plan.left),w.scratch)
    A.execute_prepared!(w,input.positions,input.cell,ps)
    owning=A.prepared_result(w); saved=deepcopy(owning)
    for algorithm in (:optimized,:dag,:reference,:dag)
        A.execute_prepared!(w,input.positions,input.cell,ps;product_algorithm=algorithm)
        @test A.prepared_result(w).forces ≈ saved.forces atol=2e-10 rtol=2e-10
    end
    @test owning==saved
    @test w.prepared.product_plan===plan
    legacy=A.PreparedACE{typeof(m),typeof(p.basis),typeof(p.weights),4}(m,p.basis,p.weights,true,false)
    @test legacy.product_algorithm===:auto
    legacy6=typeof(p)(m,p.basis,p.weights,true,false,:dag)
    @test legacy6.product_plan.terminals==plan.terminals
    legacy_s=A.SiteScratch((getfield(w.scratch[1],k) for k in 1:15)...)
    site=A.NL.SiteView(input.environment,input.positions,input.cell,1)
    A.NI.evaluate_site!(w.slots[1].outputs[1],p,site,ps,legacy_s,(A.AC.Energy(),A.AC.Forces(),A.AC.Virial()))
    @test length(legacy_s.dag_values)==length(w.scratch[1].dag_values)
    generic=A.NI.prepare_potential(m,[AtomsBase.species(sys,i) for i=1:length(sys)];width=4,product_algorithm=:dag)
    A.NI.update_prepared!(generic,ps)
    @test generic.product_algorithm===:dag
    @test length(A.NI.make_scratch(generic,0,Float32).dag_values)==length(w.scratch[1].dag_values)
    # A changed structural model gets a new plan; old prepared metadata rejects
    # in-place mutation before execution rather than silently using stale nodes.
    old=copy(m.products[1]); m.products[1][1]=mod1(old[1]+1,length(m.density))
    @test_throws ArgumentError A.execute_prepared!(w,input.positions,input.cell,ps)
    @test_throws ArgumentError A.prepare_executor!(w,input.environment)
    @test_throws ArgumentError A.NI.update_prepared!(generic,ps)
    m.products[1].=old
end

@testset "DAG prepared integration: pruning, custom radials, growth and batching" begin
    A=ACECore; NL=A.NL
    m,ps,_=A.load_resolved(ACE_MODEL_PATH)
    selection=A.BasisSelection(m;indices=[1,4,9,17,33,57])
    pruned=A.prune_basis(m,ps,selection)
    empty=A.prune_basis(m,ps,A.BasisSelection(m;indices=Int[]))
    metadata=A.RadialMetadata(m.radial_spec;cutoff=5.5,simd=true,ad=true,origin=:regular)
    ch=A.CallableRadial(A.SmoothExponential(),metadata;config=(powers=[n-1 for (n,l) in m.radial_spec],cutoff=5.5))
    custom=A.replace_radial(m,ps,fill(ch,2,2);parameters=(rate=[.7],amplitude=[.2]))
    systems=(TestSystem(zeros(3,1),9Matrix{Float64}(I,3,3),(false,false,false),[:Ti]),
        TestSystem([.1 2.3 1.1;.2 .7 2.4;.3 1.2 1.5],[3.7 .7 .1;0 3.8 .3;0 0 3.6],(true,true,true),[:Ti,:Al,:Ti]),
        TestSystem([0. 2.5 15. 17.5 30. 32.5;zeros(2,6)],50Matrix{Float64}(I,3,3),(false,false,false),[:Ti,:Al,:Ti,:Al,:Ti,:Al]))
    for bundle in ((model=m,parameters=ps),pruned,empty,custom), T in (Float32,Float64), width in (1,4)
        b=A.adapt_model(bundle.model,bundle.parameters;precision=T)
        env=A.FullEnvironment(;width); w=A.CPUExecutor(b.model,env;precision=T,slots=2,product_algorithm=:dag)
        graph=w.prepared.product_plan
        for sys in (systems...,systems[1],systems[2])
            input=A.ace_inputs(sys,b.model;environment=env,precision=T)
            A.prepare_executor!(w,env)
            @test w.prepared.product_plan===graph
            for simd in (false,true), centers in (nothing,[1])
                expected=A.ace_evaluate(input.positions,b.model,b.parameters,env,input.cell;centers)
                A.execute_prepared!(w,input.positions,input.cell,b.parameters;simd,centers,parallel=true)
                y=A.prepared_result(w); tol=T===Float32 ? 3e-4 : 2e-10
                @test y.energy ≈ expected.energy atol=tol rtol=tol
                @test y.forces ≈ expected.forces atol=tol rtol=tol
                @test y.virial ≈ expected.virial atol=tol rtol=tol
            end
        end
        for batch_size in (1,2,4)
            generic=NL.prepare_site_execution(b.model,systems[3];precision=T,width,slots=2,batch_size,product_algorithm=:dag)
            NL.execute!(generic.execution,generic.positions,generic.cell,b.parameters;parallel=true)
            y=NL.site_result(generic.execution)
            expected=A.ace_evaluate(generic.positions,b.model,b.parameters,generic.execution.environment,generic.cell)
            tol=T===Float32 ? 3e-4 : 2e-10
            @test y.energy ≈ expected.energy atol=tol rtol=tol
            @test y.forces ≈ expected.forces atol=tol rtol=tol
            @test y.virial ≈ expected.virial atol=tol rtol=tol
        end
    end
end
