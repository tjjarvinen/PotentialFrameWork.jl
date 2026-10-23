using SpheriCart, TOML, LinearAlgebra
const ACECore=PotentialFrameWork.ACE
const ACE_MODEL_PATH=joinpath(@__DIR__,"../data/ace_ref/resolved_model_v1.toml")

@testset "resolved ACE schema and parameters" begin
    m,ps,st=ACECore.load_resolved(ACE_MODEL_PATH)
    @test m.atomic_numbers==[22,13]
    @test size(ps.WB)==(123,2)
    @test size(ps.Wpair)==(12,2)
    @test length(m.density)==27
    @test size(m.coupling)==(123,167)
    @test st==NamedTuple()
    for invalid in (Dict("schema_version"=>2),Dict("schema"=>"PotentialFrameWork.resolved_ACE","schema_version"=>2))
        mktemp() do path,io
            TOML.print(io,invalid); close(io)
            @test_throws ArgumentError ACECore.load_resolved(path)
        end
    end
end

@testset "complete species environments, padding and write footprints" begin
    model,ps,st=ACECore.load_resolved(ACE_MODEL_PATH)
    H=[2.1 0.6 0.2;0 2.2 0.4;0 0 2.0]
    for width in (1,2,4,8,16), pb in ((true,true,true),(true,false,true),(false,false,false))
        env=ACECore.FullEnvironment(;width)
        for n in (3,1,0,2)
            xyz=H*[0.1 0.4 0.8;0.2 0.6 0.3;0.1 0.7 0.5][:,1:n]
            n>0 && (xyz[:,1]+=2H[:,1])
            s=TestSystem(xyz,H,pb,[:Ti,:Al,:Ti][1:n])
            input=ACECore.ace_inputs(s,model;environment=env)
            for i=1:n
                actual=Set{Tuple{Int,NTuple{3,Int}}}()
                for z=1:2
                    range=env.offsets[z,i]:env.offsets[z+1,i]-1
                    @test length(range)%width==0
                    for k in range
                        d=ACECore.edge_displacement(input.positions,input.cell,env,i,k)
                        @test all(isfinite,d)
                        if env.active[k]
                            j=env.neighbours[k]
                            @test env.species[j]==z
                            @test j in env.footprints[i]
                            @test d ≈ xyz[:,j]+H*env.shifts[k]-xyz[:,i] atol=1e-12
                            push!(actual,(j,Tuple(env.shifts[k])))
                        else
                            @test d==[1,0,0]
                        end
                    end
                end
                expected=Set{Tuple{Int,NTuple{3,Int}}}()
                for j=1:n, x=-6:6,y=-6:6,z=-6:6
                    shift=(x,y,z)
                    any(a -> !pb[a] && shift[a]!=0,1:3) && continue
                    i==j && shift==(0,0,0) && continue
                    norm(xyz[:,j]+H*collect(shift)-xyz[:,i])<5.5 || continue
                    push!(expected,(j,shift))
                end
                @test actual==expected
            end
            @test sort(vcat(Int[],env.phases...))==collect(1:n)
            for phase in env.phases
                writes=Int[]
                for i in phase; append!(writes,env.footprints[i]); end
                @test length(unique(writes))==length(writes)
            end
        end
    end
    for label in (ChemicalSpecies(:Ti;atom_name=:TI),ChemicalSpecies(:Ti48),ChemicalSpecies(:H))
        s=TestSystem(zeros(3,1),H,(false,false,false),[label])
        @test_throws ArgumentError ACECore.ace_inputs(s,model)
    end
end

@testset "independent ACE derivatives and invariance" begin
    model,ps,st=ACECore.load_resolved(ACE_MODEL_PATH)
    H=[4.7 0.5 0.1;0 4.8 0.3;0 0 4.6]
    s=TestSystem([0.1 2.3 1.1;0.2 0.7 2.4;0.3 1.2 1.5],H,(true,true,true),[:Ti,:Al,:Ti])
    inp=ACECore.ace_inputs(s,model); X=inp.positions; env=inp.environment
    eval(x,h,params=ps)=ACECore.ace_evaluate(x,model,params,env,h)
    r=eval(X,H)
    @test r.forces ≈ -central_gradient(x -> eval(x,H).energy,X) atol=2e-6 rtol=1e-6
    @test sum(r.forces;dims=2) ≈ zeros(3,1) atol=1e-11
    @test eval(X .+ [0.3,-0.7,0.1],H).energy ≈ r.energy atol=1e-10
    I3=Matrix{Float64}(I,3,3); h=2e-6
    for strain in (I3,[0. 1 0;0 0 0;0 0 0.],[0. 0 0;0 0 1;0 0 0.])
        plus=I3+h*strain; minus=I3-h*strain
        fd=(eval(plus*X,plus*H).energy-eval(minus*X,minus*H).energy)/(2h)
        @test fd ≈ -sum(r.virial.*transpose(strain)) atol=2e-6 rtol=1e-6
    end
    rot=[0. -1 0;1 0 0;0 0 1.]
    rr=eval(rot*X,rot*H)
    @test rr.energy ≈ r.energy atol=1e-10
    @test rr.forces ≈ rot*r.forces atol=1e-10
    @test rr.virial ≈ rot*r.virial*transpose(rot) atol=1e-10
    zero_ps=(WB=zero(ps.WB),Wpair=zero(ps.Wpair))
    rz=eval(X,H,zero_ps)
    @test rz.energy==sum(model.reference[env.species])
    @test all(iszero,rz.forces) && all(iszero,rz.virial)
    @test ACECore.ace_evaluate(X,model,ps,env,H;energy_only=true).energy==r.energy
    design=ACECore.design_matrix(X,model,ps,env,H)
    predicted=design.matrix*ACECore.linear_parameters(ps)+design.baseline
    @test predicted ≈ vcat(r.energy,vec(r.forces),vec(r.virial)) atol=1e-9 rtol=1e-10
    fixed_pair=ACECore.design_matrix(X,model,ps,env,H;fit_pair=false)
    @test fixed_pair.matrix*vec(ps.WB)+fixed_pair.baseline ≈ predicted atol=1e-9
    @test ACECore.install_readout(ps,vec(ps.WB);fit_pair=false).Wpair==ps.Wpair
    selected=ACECore.training_rows(design,(energy=r.energy,forces=r.forces,virial=nothing);weights=(energy=3.,forces=0.4,virial=1.))
    @test selected.matrix*ACECore.linear_parameters(ps) ≈ selected.labels atol=1e-9
    @test selected.row_weights==vcat(3.,fill(0.4,length(X)))
    @test size(ACECore.smoothness_prior(model))==(270,270)
    @test all(>(0),diag(ACECore.smoothness_prior(model)))
    coeff=0.3ACECore.linear_parameters(ps)
    installed=ACECore.install_readout(ps,coeff)
    @test ACECore.linear_parameters(installed)==coeff
    @test ACECore.linear_parameters(ps)!=coeff
    ri=eval(X,H,installed)
    @test design.matrix*coeff+design.baseline ≈ vcat(ri.energy,vec(ri.forces),vec(ri.virial)) atol=1e-9
    calc=ACECore.ACECalculator(model,ps)
    AC.Testing.test_energy_forces_virial(s,calc)
    expected=AC.energy_forces_virial(s,calc)
    @test ustrip(expected.energy) ≈ r.energy
    @test ustrip.(force_matrix(expected.forces)) ≈ r.forces
    @test ustrip.(expected.virial) ≈ r.virial
    f=AC.zero_forces(s,calc)
    AC.energy_forces_virial!(f,s,calc)
    AC.forces!(f,s,calc)
    @test force_matrix(f) ≈ 2force_matrix(expected.forces)
    AC.set_parameters!(calc,installed)
    @test ustrip(AC.potential_energy(s,calc)) ≈ ri.energy
    @test AC.calculate((AC.Energy(),AC.Forces(),AC.Virial()),s,calc;parameters=ps,state=(token=3,)).state==(token=3,)
    empty=ACECore.ace_inputs(TestSystem(zeros(3,0),H,s.pbc,Symbol[]),model)
    re=ACECore.ace_evaluate(empty.positions,model,ps,empty.environment,empty.cell)
    @test re.energy==0 && size(re.forces)==(3,0) && all(iszero,re.virial)
end
