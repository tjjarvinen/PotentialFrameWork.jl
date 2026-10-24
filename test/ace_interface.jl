using TOML, SparseArrays
@testset "ACE interface: selection, centers and species orders" begin
    A=ACECore
    m,readout,st=A.load_resolved(ACE_MODEL_PATH)
    ps=A.trainable_parameters(m,readout)
    ps.radial.scale .*= 0.83; ps.pair.scale .*= 1.07
    H=[4.7 0.5 0.1;0 4.8 0.3;0 0 4.6]
    sys=TestSystem([0.1 2.3 1.1;0.2 0.7 2.4;0.3 1.2 1.5],H,(true,true,true),[:Ti,:Al,:Ti])
    inp=A.ace_inputs(sys,m;environment=A.FullEnvironment(width=4)); X=inp.positions; env=inp.environment
    evaluate(model,params=ps;kwargs...)=A.ace_evaluate(X,model,params,env,H;kwargs...)
    function same(a,b;atol=1e-9)
        @test a.energy ≈ b.energy atol=atol
        @test a.forces ≈ b.forces atol=atol
        @test a.virial ≈ b.virial atol=atol
    end
    original=evaluate(m)
    typed=A.ResolvedModel{Float64}(m.atomic_numbers,m.reference,m.radial,m.pair_radial,m.radial_spec,m.pair_spec,m.lmax,m.density,m.products,m.coupling,m.source_sha256)
    same(evaluate(typed),original)
    ds=A.descriptors(m)
    @test length(ds)==123 && ds[1].original_index==1
    @test ds[1].nominal_body_orders==ds[1].orders .+ 1
    @test !isnothing(ds[1].index_weight)
    @test A.basis_summary(m).dimensions.products==167
    s=A.BasisSelection(m;indices=[1,4,9,17,33,57])
    masked=A.mask_basis(m,s); pruned=A.prune_basis(m,ps,s)
    @test length(pruned.model.products)<length(m.products)
    @test length(pruned.model.density)<length(m.density)
    @test length(pruned.model.radial_spec)<length(m.radial_spec)
    @test A.descriptors(pruned.model)[3].id==ds[9].id
    twice=A.prune_basis(pruned.model,pruned.parameters,A.BasisSelection(pruned.model;indices=[3,1]))
    @test A.descriptors(twice.model)[1].id==ds[9].id
    @test [t.original_density for t in A.descriptors(twice.model)[1].terms]==[t.original_density for t in ds[9].terms]
    @test [t.original_radial for t in A.descriptors(twice.model)[1].terms]==[t.original_radial for t in ds[9].terms]
    @test_throws ArgumentError A.mask_basis(pruned.model,s)
    @test_throws ArgumentError A.BasisSelection(m;indices=[1,1])
    same(evaluate(masked),evaluate(pruned.model,pruned.parameters))
    mapped_gradient=Zygote.gradient(p -> A.ace_evaluate(X,pruned.model,A.restrict_parameters(p,pruned.maps),env,H).energy,ps)[1]
    masked_gradient=Zygote.gradient(p -> A.ace_evaluate(X,masked,p,env,H).energy,ps)[1]
    @test PP.parameter_flatten(mapped_gradient) ≈ PP.parameter_flatten(masked_gradient) atol=1e-10
    same(evaluate(A.mask_basis(masked,A.BasisSelection(masked))),original)
    coeff=PP.parameter_flatten(ps); direction=sin.(1:length(coeff))
    updated=PP.parameter_rebuild(ps,coeff+0.1direction)
    reduced=A.prune_basis(m,updated,s)
    same(evaluate(masked,updated),evaluate(reduced.model,reduced.parameters))
    @test A.linear_parameters(pruned.parameters) ≈ transpose(pruned.coefficient_map)*A.linear_parameters(ps)
    R=A.smoothness_prior(m)
    @test A.mapped_prior(R,pruned.coefficient_map)*A.linear_parameters(pruned.parameters) ≈ R*(pruned.coefficient_map*A.linear_parameters(pruned.parameters))
    for centers in (nothing,1:3,[3,1],Int[])
        result=evaluate(m;centers)
        indices=isnothing(centers) ? (1:3) : centers
        @test result.energy ≈ sum((A.site_diagnostics(X,m,ps,env,H,i).energy for i in indices);init=0.) atol=1e-10
        w=A.CPUExecutor(m,env;slots=2)
        A.execute_prepared!(w,X,H,ps;centers,parallel=true)
        same(A.prepared_result(w),result)
        @test sizeof(w.force)==3length(sys)*sizeof(Float64)
    end
    for c in ([1,1],[0],[4],[true])
        @test_throws ArgumentError evaluate(m;centers=c)
    end
    @test any(!iszero,evaluate(m;centers=[1]).forces[:,2:3])
    single=[evaluate(m;centers=[i]) for i=1:3]
    same((energy=sum(r.energy for r in single),forces=sum(r.forces for r in single),virial=sum(r.virial for r in single)),original)
    total=zero(original.energy); f=zeros(size(X)); v=zeros(3,3)
    for k=0:3
        projected=A.project_species(m,[:Ti],k:k)
        result=evaluate(projected)
        total+=result.energy; f+=result.forces; v+=result.virial
        w=A.CPUExecutor(projected,env;slots=2)
        A.execute_prepared!(w,X,H,ps;parallel=true)
        same(A.prepared_result(w),result)
    end
    same((energy=total,forces=f,virial=v),original)
    same(evaluate(A.project_species(m,[:Ti,:Al],0:3)),original)
    same(evaluate(A.project_species(m,Symbol[],0:3)),original)
    @test all(iszero,evaluate(A.project_species(m,Symbol[],1:3)).forces)
    @test evaluate(A.project_species(m,[:Ti],9:10)).energy==0
    @test_throws ArgumentError A.project_species(m,[:H],0:1)
    @test_throws ArgumentError A.project_species(m,[:Ti],-1:2)
    @test_throws ArgumentError A.project_species(m,[:Ti],[0,1])
    model=A.project_species(masked,[:Al],1:2); centers=[3,1]
    y=evaluate(model;centers)
    @test -y.forces ≈ central_gradient(x -> A.ace_evaluate(x,model,ps,env,H;centers).energy,X) atol=3e-6
    for fit_pair in (true,false)
        d=A.design_matrix(X,model,ps,env,H;centers,fit_pair)
        coefficients=fit_pair ? A.linear_parameters(ps) : vec(ps.WB)
        @test d.matrix*coefficients+d.baseline ≈ vcat(y.energy,vec(y.forces),vec(y.virial)) atol=1e-9
    end
    B=sin.(X); C=reshape(cos.(1:9),3,3)
    loss(x,p)=begin
        r=A.ace_evaluate(x,model,p,env,H;centers)
        0.3r.energy+sum(B.*r.forces)+sum(C.*r.virial)
    end
    gx,gp=Zygote.gradient(loss,X,ps)
    @test gx ≈ central_gradient(x -> loss(x,ps),X) atol=5e-6 rtol=5e-6
    @test gx ≈ ForwardDiff.gradient(x -> loss(x,ps),X) atol=1e-9
    @test all(iszero,gp.WB[setdiff(1:123,collect(s.indices)),:])
    h=2e-6
    fd=(loss(X,PP.parameter_rebuild(ps,coeff+h*direction))-loss(X,PP.parameter_rebuild(ps,coeff-h*direction)))/(2h)
    @test dot(PP.parameter_flatten(gp),direction) ≈ fd atol=4e-6 rtol=4e-6
    empty=A.prune_basis(m,ps,A.BasisSelection(m;indices=Int[]))
    @test size(empty.model.coupling)==(0,0)
    same(evaluate(A.mask_basis(m,A.BasisSelection(m;indices=Int[]))),evaluate(empty.model,empty.parameters))
end

@testset "ACE custom radials, persistence and force-loss derivatives" begin
    A=ACECore
    m,base,st=A.load_resolved(ACE_MODEL_PATH)
    metadata=A.RadialMetadata(m.radial_spec;cutoff=5.5,index_weights=sum.(m.radial_spec),simd=true,ad=true,origin=:regular)
    config=(powers=[n-1 for (n,l) in m.radial_spec],cutoff=5.5)
    ch=A.CallableRadial(A.SmoothExponential(),metadata;config,
        persistence=(route="registered",id="PotentialFrameWork.smooth_exponential.v1"))
    custom_ps=(rate=[0.7],amplitude=[0.2])
    replaced=A.replace_radial(m,base,fill(ch,2,2);parameters=custom_ps)
    model=A.project_species(A.mask_basis(replaced.model,A.BasisSelection(replaced.model;indices=[1,4,10,35,71])),[:Ti],1:2)
    ps=A.trainable_parameters(model,replaced.parameters)
    ps.radial.scale .*= 0.9
    H=9Matrix{Float64}(I,3,3)
    sys=TestSystem([0.1 2.3 1.1;0.2 0.7 2.4;0.3 1.2 1.5],H,(false,false,false),[:Ti,:Al,:Ti])
    input=A.ace_inputs(sys,model;environment=A.FullEnvironment(width=8)); X=input.positions; env=input.environment
    caps=A.radial_capabilities(ch,custom_ps)
    @test caps.cpu_scalar && caps.cpu_simd && caps.training_ad_declared && !caps.gpu_verified
    @test caps.training_ad_checked && caps.training_ad_local
    for T in (Float32,Float64), width in (1,2,4,8)
        rs=ntuple(i -> T(0.3+0.2i),width)
        r=PotentialFrameWork.ACE.Vec(rs)
        values=fill(zero(r),length(m.radial_spec)); derivs=copy(values)
        A.radial_evaluate!(values,derivs,ch,r,custom_ps)
        for lane=1:width
            v,d=A.radial_ed(ch,rs[lane],custom_ps)
            @test [x[lane] for x in values] ≈ v rtol=3e-6 atol=1e-8
            @test [x[lane] for x in derivs] ≈ d rtol=3e-6 atol=1e-8
        end
    end
    r=1.3; h=2e-5
    v,d=A.radial_ed(ch,r,custom_ps)
    @test d ≈ (A.radial_values(ch,r+h,custom_ps)-A.radial_values(ch,r-h,custom_ps))/(2h) atol=1e-10
    @test all(iszero,A.radial_values(ch,5.5,custom_ps))
    values_only=(v,r,p,c)->A.SmoothExponential()(v,nothing,r,p,c)
    fallback=A.CallableRadial(values_only,metadata;config,derivative=:ad)
    @test last(A.radial_ed(fallback,r,custom_ps)) ≈ d atol=1e-12
    alloc=A.allocating_radial((r,p,c)->A.radial_ed(ch,r,p),metadata;config)
    alloc_ad=A.allocating_radial((r,p,c)->A.radial_values(ch,r,p),metadata;config,derivative=:ad)
    @test first(A.radial_ed(alloc,r,custom_ps)) ≈ v
    @test last(A.radial_ed(alloc_ad,r,custom_ps)) ≈ d
    @test_throws ArgumentError A.RadialMetadata(m.radial_spec;cutoff=5.5,envelope=:external)
    unknown=A.CallableRadial(A.SmoothExponential(),A.RadialMetadata(m.radial_spec;cutoff=5.5,ad=true);config)
    unknown_model=A.replace_radial(m,base,fill(unknown,2,2);parameters=custom_ps).model
    @test isnothing(first(A.descriptors(unknown_model)).index_weight)
    @test_throws ArgumentError A.smoothness_prior(unknown_model)
    @test all(==(1),diag(A.smoothness_prior(unknown_model;unknown=:neutral)))
    @test_throws ArgumentError A.BasisSelection(unknown_model;max_index_weight=3)
    unknown_pair=A.CallableRadial(A.SmoothExponential(),A.RadialMetadata(m.pair_spec;cutoff=5.5,ad=true);
        config=(powers=[n-1 for (n,l) in m.pair_spec],cutoff=5.5))
    pair_model=A.replace_radial(m,base,fill(unknown_pair,2,2);parameters=custom_ps,pair=true).model
    @test_throws ArgumentError A.smoothness_prior(pair_model)
    @test all(==(1),diag(A.smoothness_prior(pair_model;unknown=:neutral)))
    @test size(A.smoothness_prior(pair_model;fit_pair=false))==(length(base.WB),length(base.WB))
    @test_throws DimensionMismatch A.radial_evaluate!(zeros(2),zeros(2),ch,1.,custom_ps)
    @test_throws ArgumentError A.replace_radial(m,base,fill(A.CallableRadial(A.SmoothExponential(),A.RadialMetadata([(1,0)];cutoff=5.5);config),2,2);parameters=custom_ps)
    centers=[1,3]
    y=A.ace_evaluate(X,model,ps,env,H;centers)
    w=A.CPUExecutor(model,env;slots=2)
    A.execute_prepared!(w,X,H,ps;parallel=true,centers)
    @test A.prepared_result(w).forces ≈ y.forces atol=1e-11
    @test A.prepared_result(w).energy ≈ y.energy atol=1e-11
    B=sin.(X); C=reshape(cos.(1:9),3,3)
    loss(x,p)=begin
        result=A.ace_evaluate(x,model,p,env,H;centers)
        result.energy+sum(B.*result.forces)+sum(C.*result.virial)
    end
    gx,gp=Zygote.gradient(loss,X,ps)
    @test gx ≈ central_gradient(x -> loss(x,ps),X) atol=5e-6
    flat=PP.parameter_flatten(ps); direction=cos.(1:length(flat)); h=2e-6
    @test dot(PP.parameter_flatten(gp),direction) ≈ (loss(X,PP.parameter_rebuild(ps,flat+h*direction))-loss(X,PP.parameter_rebuild(ps,flat-h*direction)))/(2h) atol=4e-6 rtol=4e-6
    for key in (:rate,:amplitude)
        delta=zeros(length(flat)); gp_custom=getproperty(gp.custom_radial,key)[1]
        pplus=deepcopy(ps); pminus=deepcopy(ps)
        getproperty(pplus.custom_radial,key)[1]+=h; getproperty(pminus.custom_radial,key)[1]-=h
        @test gp_custom ≈ (loss(X,pplus)-loss(X,pminus))/(2h) atol=2e-6
    end
    _,pb=rrule(A.ace_evaluate,X,model,ps,env,H;centers)
    before=pb((energy=1.,forces=B,virial=C))
    centers[1]=2
    @test pb((energy=1.,forces=B,virial=C))[2] ≈ before[2]
    @test pb(ZeroTangent())[2] isa AbstractZero
    adapted=A.replace_radial(m,base,fill(fallback,2,2);parameters=custom_ps)
    fallback_loss(x,p)=sum(B.*A.ace_evaluate(x,adapted.model,p,env,H).forces)
    fgx,fgp=Zygote.gradient(fallback_loss,X,adapted.parameters)
    @test fgx ≈ central_gradient(x -> fallback_loss(x,adapted.parameters),X) atol=4e-6
    fp=deepcopy(adapted.parameters); fm=deepcopy(adapted.parameters)
    fp.custom_radial.rate[1]+=h; fm.custom_radial.rate[1]-=h
    @test fgp.custom_radial.rate[1] ≈ (fallback_loss(X,fp)-fallback_loss(X,fm))/(2h) atol=4e-6
    dsgn=A.design_matrix(X,model,ps,env,H;centers=[1,3])
    @test dsgn.matrix*A.linear_parameters(ps)+dsgn.baseline ≈ vcat(y.energy,vec(y.forces),vec(y.virial)) atol=1e-10
    rows=A.training_rows(dsgn,(energy=y.energy*u"eV",forces=y.forces*u"eV/Å",virial=y.virial*u"eV"))
    @test rows.matrix*A.linear_parameters(ps) ≈ rows.labels atol=1e-10
    reduced=A.prune_basis(model,ps,A.BasisSelection(model;indices=[1,4,10,35,71]))
    yr=A.ace_evaluate(X,reduced.model,reduced.parameters,env,H;centers=[1,3])
    @test yr.forces ≈ y.forces atol=1e-11
    table=A.tabulate_radial(ch,custom_ps;nodes=501)
    @test table.report.value_error<1e-6 && table.report.derivative_error<1e-5
    @test table.report.fixed_parameter_snapshot
    for r in (0.32,1.37,4.99)
        v,d=A.radial_ed(ch,r,custom_ps); vt,dt=A.radial_ed(table.channel,r)
        @test v ≈ vt atol=1e-6
        @test d ≈ dt atol=1e-5
    end
    @test_throws ArgumentError A.tabulate_radial(ch,custom_ps;nodes=3,atol=1e-15)
    mktempdir() do dir
        path=joinpath(dir,"custom.toml")
        calc=A.ACECalculator(reduced.model,reduced.parameters;slots=1)
        A.save_calculator(path,calc)
        restored=A.ACECalculator(path;slots=1)
        @test restored.parameters==reduced.parameters
        nested=A.prune_basis(reduced.model,reduced.parameters,A.BasisSelection(reduced.model;indices=[3,1]))
        nestedcalc=A.ACECalculator(nested.model,nested.parameters;slots=1)
        nestedpath=joinpath(dir,"nested.toml"); A.save_calculator(nestedpath,nestedcalc)
        nestedloaded=A.ACECalculator(nestedpath;slots=1)
        @test hcat(AC.forces(sys,nestedloaded;centers=[1,3])...) ≈ hcat(AC.forces(sys,nestedcalc;centers=[1,3])...) atol=1e-11u"eV/Å"
        @test first(A.descriptors(nestedloaded.model)).id==A.descriptors(reduced.model)[3].id
        rebound=A.replace_radial(reduced.model,reduced.parameters,reduced.model.radial;parameters=custom_ps)
        @test A.ace_evaluate(X,rebound.model,rebound.parameters,env,H;centers=[1,3]).energy ≈ y.energy atol=1e-11
        @test A.descriptors(restored.model)[2].id==A.descriptors(reduced.model)[2].id
        expected=AC.energy_forces_virial(sys,calc;centers=[3,1])
        actual=AC.energy_forces_virial(sys,restored;centers=[3,1])
        @test actual.energy ≈ expected.energy atol=1e-10u"eV"
        @test actual.virial ≈ expected.virial atol=1e-10u"eV"
        intact=read(path)
        unsupported=A.replace_radial(m,base,fill(unknown,2,2);parameters=custom_ps)
        @test_throws ArgumentError A.save_calculator(path,A.ACECalculator(unsupported.model,unsupported.parameters;slots=1))
        @test read(path)==intact
        tabmodel=A.replace_radial(m,base,fill(table.channel,2,2);parameters=NamedTuple())
        A.save_calculator(path,A.ACECalculator(tabmodel.model,tabmodel.parameters;slots=1))
        tabloaded=A.ACECalculator(path;slots=1)
        @test AC.potential_energy(sys,tabloaded) ≈ AC.potential_energy(sys,A.ACECalculator(tabmodel.model,tabmodel.parameters;slots=1))
        source="""
        struct SavedRadial end
        function (f::SavedRadial)(v,d,r,p,c)
            rc=c.cutoff; q=max(one(r)-r/rc,zero(r)); ex=exp(-p.rate[1]*r); q2=q*q; q3=q2*q; q4=q2*q2
            for i in eachindex(v)
                v[i]=p.amplitude[1]*q4*ex/(i+1)
                isnothing(d) || (d[i]=p.amplitude[1]*ex*(-4q3/rc-p.rate[1]*q4)/(i+1))
            end
        end
        factory(c,m)=SavedRadial()
        """
        source_ch=A.source_radial(source,"factory",metadata;config,allow_code=true)
        sr=A.replace_radial(m,base,fill(source_ch,2,2);parameters=custom_ps)
        sourcecalc=A.ACECalculator(sr.model,sr.parameters;slots=1)
        A.save_calculator(path,sourcecalc)
        @test_throws ArgumentError A.ACECalculator(path;slots=1)
        loaded=A.ACECalculator(path;slots=1,allow_code=true)
        @test AC.energy_forces_virial(sys,loaded).energy ≈ AC.energy_forces_virial(sys,sourcecalc).energy
        # Source methods were created inside this test function: public execution
        # boundaries must handle world age without invokelatest in radial loops.
        inp=A.ace_inputs(sys,loaded.model)
        @test A.ace_evaluate(inp.positions,loaded.model,loaded.parameters,inp.environment,inp.cell).energy ≈ ustrip(u"eV",AC.potential_energy(sys,loaded))
        source_loss(x,p)=sum(B.*A.ace_evaluate(x,loaded.model,p,inp.environment,inp.cell).forces)
        source_gx,source_gp=Zygote.gradient(source_loss,inp.positions,loaded.parameters)
        @test source_gx ≈ central_gradient(x -> source_loss(x,loaded.parameters),inp.positions) atol=3e-6
        sp=deepcopy(loaded.parameters); sm=deepcopy(loaded.parameters)
        sp.custom_radial.rate[1]+=h; sm.custom_radial.rate[1]-=h
        @test source_gp.custom_radial.rate[1] ≈ (source_loss(inp.positions,sp)-source_loss(inp.positions,sm))/(2h) atol=3e-6
        # Fresh-process, relocated source checkpoint: no definitions from this
        # session are available to the child process.
        moved=joinpath(dir,"relocated"); mkpath(moved)
        freshpath=joinpath(moved,"source.toml"); cp(path,freshpath)
        script="""
        using PotentialFrameWork, AtomsBase, Unitful, SpheriCart
        import AtomsCalculators as AC
        sys=AtomsBase.isolated_system([:Ti=>[0.1,0.2,0.3]*u"Å",:Al=>[2.3,0.7,1.2]*u"Å",:Ti=>[1.1,2.4,1.5]*u"Å"])
        calc=PotentialFrameWork.ACE.ACECalculator(ARGS[1];allow_code=true,slots=1)
        print(ustrip(u"eV",AC.potential_energy(sys,calc)))
        """
        fresh=parse(Float64,read(`$(Base.julia_cmd()) --startup-file=no --threads=1 --project=$(dirname(Base.active_project())) -e $script $freshpath`,String))
        @test fresh ≈ ustrip(u"eV",AC.potential_energy(sys,loaded))
        @test !isdefined(Main,:SavedRadial)
        sentinel=joinpath(dir,"sentinel")
        saved=TOML.parsefile(path)
        with_effect="write("*repr(sentinel)*", \"executed\")\n"*source
        for channel in saved["model"]["radial"]["channels"]
            channel["source"]=with_effect
            channel["source_sha256"]=bytes2hex(A.SHA.sha256(with_effect))
        end
        open(io -> TOML.print(io,saved),path,"w")
        @test_throws ArgumentError A.ACECalculator(path;slots=1)
        @test !isfile(sentinel)
        A.ACECalculator(path;allow_code=true,slots=1)
        @test read(sentinel,String)=="executed"
        rm(sentinel)
        saved["model"]["radial"]["channels"][1]["dependencies"]=Dict("MissingACERadialDependency999"=>"1.0.0")
        open(io -> TOML.print(io,saved),path,"w")
        @test_throws ArgumentError A.ACECalculator(path;allow_code=true,slots=1)
        @test !isfile(sentinel)
        # Restore a code bundle without side effects for corruption checks.
        A.save_calculator(path,sourcecalc)
        for field in ("entrypoint","interface_version","reference")
            bad=TOML.parsefile(path)
            channel=bad["model"]["radial"]["channels"][1]
            if field=="entrypoint"; channel[field]="missing_factory"
            elseif field=="interface_version"; channel[field]=999
            else; channel[field][1]["values"][1]+=1.
            end
            badpath=joinpath(dir,"bad.toml"); open(io -> TOML.print(io,bad),badpath,"w")
            @test_throws ArgumentError A.ACECalculator(badpath;allow_code=true,slots=1)
        end

        @test_throws ArgumentError A.source_radial("error(\"must not execute\")","factory",metadata)
        corrupt=TOML.parsefile(path)
        corrupt["model"]["radial"]["channels"][1]["source_sha256"]="invalid"
        open(io -> TOML.print(io,corrupt),path,"w")
        @test_throws ArgumentError A.ACECalculator(path;allow_code=true)
    end
end

@testset "independent species polynomial oracle and symmetry" begin
    A=ACECore
    products=[[1,2],[1,1,2],[2],[1,3,1]]
    for U in ([0.,2.,-0.3],[1.2,-0.7,0.2]), V in ([3.,0.,0.4],[0.,-0.2,1.1]), k=0:4
        # Independently enumerate the 2^nu factor assignments in these tiny
        # oracle products; production uses truncated convolution instead.
        expected=[sum((prod(((bits >> (t-1)) & 1)==1 ? V[a] : U[a] for (t,a) in enumerate(indices))
            for bits=0:2^length(indices)-1 if count_ones(bits)==k);init=0.) for indices in products]
        P=zeros(4); dU=zeros(3); dV=zeros(3); fwd=zeros(k+1,4); rev=similar(fwd)
        A.product_projection!(P,dU,dV,U,V,products,ones(4),k:k,fwd,rev)
        @test P ≈ expected atol=1e-13
    end
    m,ps,st=A.load_resolved(ACE_MODEL_PATH)
    H=[3.7 0.6 0.1;0 3.8 0.2;0 0 3.6]
    sys=TestSystem([0.1 2.3 1.1;0.2 0.7 2.4;0.3 1.2 1.5],H,(true,true,true),[:Ti,:Al,:Ti])
    input=A.ace_inputs(sys,m); X=input.positions; env=input.environment
    exclude=deepcopy(env)
    for i in eachindex(env.species),k=env.offsets[1,i]:env.offsets[end,i]-1
        env.species[env.neighbours[k]]==1 && (exclude.active[k]=false)
    end
    zero_order=A.ace_evaluate(X,A.project_species(m,[:Ti],0:0),ps,env,H)
    excluded=A.ace_evaluate(X,m,ps,exclude,H)
    @test zero_order.energy ≈ excluded.energy atol=1e-10
    @test zero_order.forces ≈ excluded.forces atol=1e-10
    @test zero_order.virial ≈ excluded.virial atol=1e-10
    full=A.ace_evaluate(X,m,ps,env,H)
    presence=A.ace_evaluate(X,A.project_species(m,[:Ti],1:3),ps,env,H)
    @test presence.energy ≈ full.energy-excluded.energy atol=1e-10
    @test presence.forces ≈ full.forces-excluded.forces atol=1e-10
    @test presence.virial ≈ full.virial-excluded.virial atol=1e-10
    proj=A.project_species(m,[:Ti,:Al],1:2)
    r=A.ace_evaluate(X,proj,ps,env,H;centers=[1,3])
    rotation=[0. -1 0;1 0 0;0 0 1.]
    rotated=A.ace_evaluate(rotation*X,proj,ps,env,rotation*H;centers=[1,3])
    @test rotated.energy ≈ r.energy atol=1e-10
    @test rotated.forces ≈ rotation*r.forces atol=1e-10
    @test rotated.virial ≈ rotation*r.virial*transpose(rotation) atol=1e-10
    @test sum(r.forces;dims=2) ≈ zeros(3,1) atol=1e-11
    strain=[0. 1 0;0 0 0;0 0 0.]; identity=Matrix{Float64}(I,3,3); h=2e-6
    plus=identity+h*strain; minus=identity-h*strain
    fd=(A.ace_evaluate(plus*X,proj,ps,env,plus*H;centers=[1,3]).energy-A.ace_evaluate(minus*X,proj,ps,env,minus*H;centers=[1,3]).energy)/(2h)
    @test fd ≈ -sum(r.virial.*transpose(strain)) atol=3e-6
    # A real sparse coupling row can combine different correlation orders.
    lengths=length.(m.products); j1=findfirst(==(1),lengths); j2=findfirst(==(2),lengths)
    mixed=A.ResolvedModel(m.atomic_numbers,m.reference,m.radial,m.pair_radial,m.radial_spec,m.pair_spec,m.lmax,m.density,
        m.products,sparse([1,1],[j1,j2],[1.,0.3],1,length(m.products)),m.source_sha256)
    @test only(A.descriptors(mixed)).orders==[1,2]
    @test isempty(A.BasisSelection(mixed;orders=1:1).indices)
    @test A.BasisSelection(mixed;orders=1:2).indices==(1,)
    restricted=A.restrict_basis(ACE_MODEL_PATH;indices=[9,4,1])
    reference=A.prune_basis(m,ps,A.BasisSelection(m;indices=[9,4,1]))
    @test restricted.model.coupling==reference.model.coupling
    @test restricted.parameters==reference.parameters
end
