using StaticArrays, LinearAlgebra
import AtomsCalculators as AC
const ChunkNL=PotentialFrameWork.NeighbourList
const ChunkNI=PotentialFrameWork.NeighbourInterface
const ChunkACE=PotentialFrameWork.ACE
chunk_numeric(x::AbstractArray)=map(chunk_numeric,x)
chunk_numeric(x::Number)=ustrip(x)

@testset "serial chunks: destination mapping, capacity and transpose" begin
    # Centers 1,2 occupy one cell. Their three-destination footprints have a
    # four-destination union, forcing a split at capacity three.
    sys=TestSystem([1.4 1.5 0.5 2.4;0 0 0 0;0 0 0 0],10.0Matrix(I,3,3),(true,false,false),fill(:H,4))
    n=ChunkNL.prepare_neighbors(sys,1.0;width=4)
    w=ChunkNL.form_numa_workers(n,1;scatter_capacity=3)
    @test length(w.chunks)>length(n.ranges)
    @test n.p!=collect(1:4)
    @test_throws ArgumentError ChunkNL.form_numa_workers(n,1;scatter_capacity=2)
    @test_throws ArgumentError ChunkNL.form_numa_workers(n,2)
    @test_throws ArgumentError ChunkNL.foreach_worker!(identity,w,2)
    gbar=[SVector(sin(i),cos(i),0.2i) for i=1:4]
    sorted_bar=gbar[n.p]
    lhs=Ref(0.0); rhs=Ref(0.0); site_lhs=Ref(0.0); site_rhs=Ref(0.0)
    scratch_calls=Ref(0); scratch_id=Ref{Any}()
    for repeat=1:2
        lhs[]=rhs[]=site_lhs[]=site_rhs[]=0.0
        ChunkNL.foreach_worker!(w,1;scatter=(forces=SVector{3,Float64},position_bar=SVector{3,Float64}),
            sum_outputs=(energy=Float64,vector=SVector{2,Float64}),per_site=(value=Float64,),
            init_scratch=()->(scratch_calls[]+=1; zeros(3))) do worker
            buffers=ChunkNL.get_buffers(worker)
            if !isassigned(scratch_id); scratch_id[]=ChunkNL.get_scratch(worker); end
            @test ChunkNL.get_scratch(worker)===scratch_id[]
            @test all(iszero,buffers.scatter.forces)
            @test iszero(buffers.sum_outputs.energy[])
            @test keys(buffers.per_site)==(:value,)
            for site in worker
                i=ChunkNI.central_index(site); value=0.3i
                buffers.per_site.value[ChunkNL.site_ordinal(site)]=value
                site_lhs[]+=value*(i+0.2)
                site_rhs[]+=buffers.per_site.value[ChunkNL.site_ordinal(site)]*(i+0.2)
                buffers.sum_outputs.energy[]+=1
                buffers.sum_outputs.vector[]+=SVector(1.,2.)
                for edge=1:ChunkNI.neighbour_count(site)
                    if !ChunkNI.active(site,edge)
                        @test ChunkNL.neighbor_slot(site,edge)==0
                        continue
                    end
                    d=ChunkNI.displacement(site,edge)
                    j=site.site.environment.neighbours[ChunkNL.edge_index(site.site,edge)]
                    g=0.7d
                    buffers.scatter.forces[ChunkNL.central_slot(site)]+=g
                    buffers.scatter.forces[ChunkNL.neighbor_slot(site,edge)]-=g
                    buffers.scatter.position_bar[ChunkNL.central_slot(site)]+=2g
                    buffers.scatter.position_bar[ChunkNL.neighbor_slot(site,edge)]-=2g
                    center=w.neighbours.inverse[i]
                    neighbor=worker.chunk.destinations[ChunkNL.neighbor_slot(site,edge)]
                    lhs[]+=dot(g,sorted_bar[center]-sorted_bar[neighbor])
                end
            end
        end
        r=ChunkNL.get_results(w)
        @test r.energy==4
        @test r.vector==SVector(4.,8.) # sum vectors are not atom fields
        @test r.value≈0.3collect(1:4)
        @test sum(dot.(r.forces,gbar))≈lhs[] atol=1e-14
        @test site_lhs[]≈sum(r.value .* (collect(1:4).+0.2))
        @test r.position_bar≈2 .* r.forces
        @test w.observed_capacity<=3
        scratch_id=Ref{Any}()
    end
    @test scratch_calls[]==2
    @test_throws ArgumentError ChunkNL.foreach_worker!(identity,w,1;scatter=(x=Float64,),per_site=(x=Float64,))
    @test_throws ArgumentError ChunkNL.foreach_worker!(identity,w,1;site_inputs=(:missing,))
end

function first_feature_pass!(w)
    ChunkNL.foreach_worker!(w,1;per_site=(layer1=SVector{2,Float64},)) do worker
        b=ChunkNL.get_buffers(worker)
        @test isempty(b.scatter) && isempty(b.sum_outputs)
        for site in worker
            i=ChunkNI.central_index(site)
            b.per_site.layer1[ChunkNL.site_ordinal(site)]=SVector(Float64(i),0.1i^2)
        end
    end
end
@testset "two feature passes, image aliases, publication and lifetime" begin
    sys=TestSystem([0.1 1.7 0.8;0.2 0.3 1.8;0 0 0],2.0Matrix(I,3,3),(true,true,false),[:H,:He,:H])
    n=ChunkNL.prepare_neighbors(sys,2.4;width=4)
    w=ChunkNL.form_numa_workers(n,1;scatter_capacity=3)
    @test any(i -> length(unique(n.environment.neighbours[n.environment.offsets[1,i]:n.environment.offsets[end,i]-1]))<count(n.environment.active[n.environment.offsets[1,i]:n.environment.offsets[end,i]-1]),1:3)
    first_feature_pass!(w)
    layer1=ChunkNL.export_field(w,:layer1)
    @test layer1==[SVector(Float64(i),0.1i^2) for i=1:3]
    @test w.fields[:layer1]==layer1[n.p]
    # Independent all-image reference, including repeated/self images.
    expected=copy(layer1)
    for i=1:3,j=1:3,a=-2:2,b=-2:2
        i==j && a==b==0 && continue
        d=sys.xyz[:,j]-sys.xyz[:,i]+sys.H*[a,b,0]
        norm(d)<2.4 || continue
        expected[i]+=layer1[j]*exp(-dot(d,d))
    end
    ChunkNL.foreach_worker!(w,1;site_inputs=(:layer1,),per_site=(layer2=SVector{2,Float64},)) do worker
        buffers=ChunkNL.get_buffers(worker)
        for site in worker
            value=ChunkNL.central_feature(site,Val(:layer1))
            for edge=1:ChunkNI.neighbour_count(site)
                ChunkNI.active(site,edge) || continue
                d=ChunkNI.displacement(site,edge)
                value+=ChunkNL.neighbor_feature(site,edge,Val(:layer1))*exp(-dot(d,d))
            end
            buffers.per_site.layer2[ChunkNL.site_ordinal(site)]=value
        end
    end
    @test ChunkNL.export_field(w,:layer2)≈expected
    @test ChunkNL.export_field(w,:layer1)==layer1
    @test_throws ArgumentError ChunkNL.foreach_worker!(identity,w,1;site_inputs=(:layer1,),per_site=(layer1=Float64,))
    @test_throws ArgumentError ChunkNL.foreach_worker!(identity,w,1;site_inputs=(:layer1,),scatter=(layer1=SVector{2,Float64},))
    @test_throws ErrorException ChunkNL.foreach_worker!(w,1;per_site=(failed=Float64,)) do worker
        error("deliberate failure")
    end
    @test !haskey(w.fields,:failed)
    @test_throws ArgumentError ChunkNL.foreach_worker!(identity,w,1;per_site=(incomplete=Float64,))
    @test !haskey(w.fields,:incomplete)
    @test_throws ArgumentError ChunkNL.get_results(w)
    @test ChunkNL.export_field(w,:layer1)==layer1
    first_feature_pass!(w)
    @test keys(ChunkNL.get_results(w))==(:layer1,)
    exported=ChunkNL.export_field(w,:layer1); exported[1]=zero(exported[1])
    @test ChunkNL.export_field(w,:layer1)==layer1
    n.environment.version+=UInt(1)
    @test_throws ArgumentError ChunkNL.export_field(w,:layer1)
    @test isempty(w.fields)
    empty=ChunkNL.form_numa_workers(ChunkNL.prepare_neighbors(TestSystem(zeros(3,0),10.0Matrix(I,3,3),(false,false,false),Symbol[]),1.),1)
    ChunkNL.foreach_worker!(identity,empty,1;sum_outputs=(energy=Float64,),per_site=(x=Float64,))
    @test ChunkNL.get_results(empty)==(energy=0.,x=Float64[])
end

@testset "ACE chunk E/F/V, original order, explicit state isolation" begin
    model,ps,st=ChunkACE.load_resolved(joinpath(@__DIR__,"../data/ace_ref/resolved_model_v1.toml"))
    sys=TestSystem([6.2 0.4 2.8 7.5;0.3 0.6 1.2 2.5;0.4 1.1 0.2 0.7],18.0Matrix(I,3,3),(true,false,true),[:Al,:Ti,:Al,:Ti])
    for T in (Float32,Float64), width in (1,4), batch_size in (1,3)
        c=ChunkACE.ACECalculator(model,ps,st;precision=T,width,slots=1,angular_backend=:block,batch_size)
        state=ChunkACE.ACEEvaluationState(c;scatter_capacity=4)
        tol=T===Float32 ? 3e-4 : 2e-11
        for centers in (nothing,[4,1],Int[])
            actual=AC.calculate((AC.Energy(),AC.Forces(),AC.Virial()),sys,c,c.parameters,state;centers)
            ref=ChunkACE.calculator_result(sys,c;executor=:colored,centers)
            @test actual.state===state
            @test actual.energy≈ref.energy atol=tol*u"eV" rtol=tol
            @test chunk_numeric(actual.forces)≈chunk_numeric(ref.forces) atol=tol rtol=tol
            @test actual.virial≈ref.virial atol=tol*u"eV" rtol=tol
            @test state.chunks.observed_capacity<=4
            @test state.chunks.neighbours.p!=collect(1:length(sys))
        end
        for props in ((AC.Energy(),),(AC.Forces(),),(AC.Virial(),),(AC.Energy(),AC.Virial()),(AC.Forces(),AC.Energy()),(AC.Forces(),AC.Virial()))
            actual=AC.calculate(props,sys,c,c.parameters,state)
            ref=AC.calculate(props,sys,c;executor=:colored)
            for name in keys(actual)
                name===:state && continue
                @test chunk_numeric(actual[name])≈chunk_numeric(ref[name]) rtol=tol atol=tol
            end
        end
        ref=AC.energy_forces_virial(sys,c;executor=:colored)
        r=AC.energy_forces_virial(sys,c)
        AC.energy_forces_virial(TestSystem(sys.xyz[:,1:2],sys.H,sys.pbc,sys.labels[1:2]),c)
        @test chunk_numeric(r.forces)≈chunk_numeric(ref.forces) atol=tol rtol=tol
    end
    c=ChunkACE.ACECalculator(model,ps;slots=1,angular_backend=:block)
    a=ChunkACE.ACEEvaluationState(c); b=ChunkACE.ACEEvaluationState(c)
    other=TestSystem(sys.xyz .+ [0.0 0.1 0.2 0.3;0 0 0 0;0 0 0 0],sys.H,sys.pbc,sys.labels)
    ra=AC.energy_forces_virial(sys,c;executor=:colored)
    rb=AC.energy_forces_virial(other,c;executor=:colored)
    legacy=AC.energy_forces_virial(sys,c;parallel=true)
    @test chunk_numeric(legacy.forces)≈chunk_numeric(ra.forces)
    typed=typeof(c)(c.model,c.parameters,c.state,c.environment,c.executor)
    @test AC.energy_forces_virial(sys,typed).energy≈ra.energy
    # Barrier forces overlap before inference. Distinct states retain all mutable
    # geometry/coefficients/scratch; each call executes its chunks serially.
    ready=Channel{Nothing}(2); go=Channel{Nothing}(2)
    function call(system,state)
        put!(ready,nothing); take!(go)
        AC.calculate((AC.Energy(),AC.Forces(),AC.Virial()),system,c,ps,state)
    end
    ta=Threads.@spawn call(sys,a); tb=Threads.@spawn call(other,b)
    take!(ready); take!(ready); put!(go,nothing); put!(go,nothing)
    ya=fetch(ta); yb=fetch(tb)
    forwarded=AC.calculate((AC.Energy(),AC.Forces(),AC.Virial()),sys,c,ps,a;parallel=true)
    @test chunk_numeric(forwarded.forces)≈chunk_numeric(ra.forces)
    @test chunk_numeric(ya.forces)≈chunk_numeric(ra.forces)
    @test chunk_numeric(yb.forces)≈chunk_numeric(rb.forces)
    @test ya.energy≈ra.energy && yb.energy≈rb.energy
    @test ya.virial≈ra.virial && yb.virial≈rb.virial
    @test a.scratch.numerical!==b.scratch.numerical
    @test a.environment!==b.environment
    changed=merge(ps,(WB=ps.WB.*1.1,))
    y=AC.calculate((AC.Energy(),AC.Forces()),sys,c,changed,a)
    ref=AC.calculate((AC.Energy(),AC.Forces()),sys,c,changed;executor=:colored)
    @test y.energy≈ref.energy && chunk_numeric(y.forces)≈chunk_numeric(ref.forces)
    # Caller contract with the convenient calculator scratch adapter. Species
    # order is deliberately Al,Ti here, opposite to the model's channel order.
    workgroups=ChunkNL.form_numa_workers(ChunkNL.prepare_neighbors(sys,c),1)
    ChunkNL.foreach_worker!(workgroups,1;
        scatter=(forces=SVector{3,Float64},),
        sum_outputs=(energy=Float64,virial=SMatrix{3,3,Float64,9}),
        init_scratch=()->ChunkACE.create_scratch(c)) do worker
        for site in worker
            PotentialFrameWork.evaluate!(ChunkNL.get_buffers(worker),ChunkNL.get_scratch(worker),c,site)
        end
    end
    raw=ChunkNL.get_results(workgroups)
    @test raw.energy≈ustrip(ra.energy)
    @test raw.forces≈chunk_numeric(ra.forces)
    @test raw.virial≈ustrip.(ra.virial)
    # Small periodic self-image fixture retains virial despite zero atom force.
    tiny=TestSystem(reshape([0.2,0.3,0.4],3,1),2.0Matrix(I,3,3),(true,true,true),[:Ti])
    x=AC.energy_forces_virial(tiny,c)
    ref=AC.energy_forces_virial(tiny,c;executor=:colored)
    @test x.energy≈ref.energy && x.virial≈ref.virial
    @test maximum(norm.(x.forces))<1e-10u"eV/Å"
    @test !iszero(norm(x.virial))
end

@testset "chunk ACE versus saved 33-frame oracle" begin
    fixture_dir=joinpath(@__DIR__,"../data/ace_ref")
    reference=ChunkACE.TOML.parsefile(joinpath(fixture_dir,"reloaded_predictions.toml"))
    model,ps,st=ChunkACE.load_resolved(joinpath(fixture_dir,"resolved_model_v1.toml"))
    frames=TestSystem[]
    # This fixture-only reader takes geometry directly from the saved XYZ.
    open(joinpath(fixture_dir,"reference_results.xyz")) do io
        while !eof(io)
            n=parse(Int,readline(io)); header=readline(io)
            lattice=match(r"Lattice=\"([^\"]+)\"",header)
            H=reshape(parse.(Float64,split(lattice[1])),3,3)
            X=zeros(3,n); labels=Symbol[]
            for i=1:n
                fields=split(readline(io)); push!(labels,Symbol(fields[1]))
                X[:,i]=parse.(Float64,fields[2:4])
            end
            push!(frames,TestSystem(X,H,(true,true,true),labels))
        end
    end
    @test length(frames)==length(reference["frames"])==33
    for batch_size in (1,5)
        c=ChunkACE.ACECalculator(model,ps;slots=1,batch_size,angular_backend=:block)
        state=ChunkACE.ACEEvaluationState(c;scatter_capacity=64)
        for (sys,ref) in zip(frames,reference["frames"])
            y=AC.calculate((AC.Energy(),AC.Forces(),AC.Virial()),sys,c,ps,state)
            @test ustrip(y.energy)≈ref["energy_eV"] atol=1e-8 rtol=0
            @test hcat(chunk_numeric(y.forces)...)≈hcat(ref["forces_eV_per_A"]...) atol=1e-9 rtol=1e-10
            @test ustrip.(y.virial)≈hcat(ref["virial_eV_columns"]...) atol=1e-8 rtol=1e-10
            colored=ChunkACE.calculator_result(sys,c;executor=:colored)
            @test ustrip(y.energy)≈ustrip(colored.energy) atol=1e-8 rtol=0
            @test chunk_numeric(y.forces)≈chunk_numeric(colored.forces) atol=1e-9 rtol=1e-10
            @test ustrip.(y.virial)≈ustrip.(colored.virial) atol=1e-8 rtol=1e-10
        end
    end
end
