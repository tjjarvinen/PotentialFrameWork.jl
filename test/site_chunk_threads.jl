# This file shares the fixture bindings in site_chunks.jl.
using StaticArrays, LinearAlgebra

struct RejectChunkMerge
    value::Int
end
Base.zero(::Type{RejectChunkMerge})=RejectChunkMerge(0)
Base.:+(::RejectChunkMerge,::RejectChunkMerge)=error("deliberate merge failure")

@testset "Threaded chunk routing, publication and recovery" begin
    sys=TestSystem([1.4 1.5 0.5 2.4;0 0 0 0;0 0 0 0],10.0Matrix(I,3,3),
                   (true,false,false),fill(:H,4))
    n=ChunkNL.prepare_neighbors(sys,1.0;width=4)
    w=ChunkNL.form_numa_workers(n,3;scatter_capacity=3)
    serial=ChunkNL.form_numa_workers(n,1;scatter_capacity=3)
    @test length(w.groups)==3
    @test reduce(vcat,w.groups)==collect(1:4)
    @test length(w.chunks)>length(n.ranges)
    @test n.p!=collect(1:4)

    function feature_force!(layout,workers;timings=nothing)
        ChunkNL.foreach_worker!(layout,workers;scatter=(forces=SVector{3,Float64},),timings,
                                per_site=(layer=Float64,),sum_outputs=(energy=Float64,)) do worker
            b=ChunkNL.get_buffers(worker)
            for site in worker
                i=ChunkNI.central_index(site)
                b.per_site.layer[ChunkNL.site_ordinal(site)]=Float64(i)
                # Centers 1 and 2 are separate chunks. They both reach
                # unselected atoms, including across ownership boundaries.
                i in (1,2) || continue
                b.sum_outputs.energy[]+=i
                for k in 1:ChunkNI.neighbour_count(site)
                    j=ChunkNL.neighbor_slot(site,k)
                    j==0 && continue
                    d=ChunkNI.displacement(site,k)
                    b.scatter.forces[ChunkNL.central_slot(site)]+=d
                    b.scatter.forces[j]-=d
                end
            end
        end
        ChunkNL.get_results(layout)
    end

    try
        original=feature_force!(w,2)
        @test original==feature_force!(serial,1)
        @test original.energy==3 && original.layer==[1.,2.,3.,4.]
        @test all(i -> !iszero(original.forces[i]),(3,4))
        pool=w.pool
        @test length(pool.tasks)==2 && all(task -> !istaskdone(task),pool.tasks)
        @test feature_force!(w,2)==original && w.pool===pool
        storage=pool.storage
        stats=ChunkNL.ChunkTimings()
        @test feature_force!(w,2;timings=stats)==original && pool.storage===storage
        @test length(stats.site)==2 && stats.pass>=stats.merge_barrier
        @test 1<=stats.peak_slots<=2*ChunkNL.CHUNKS_PER_WAVE && stats.capacity<=w.scatter_capacity
        @test feature_force!(w,3)==original
        @test all(istaskdone,pool.tasks)

        ChunkNL.foreach_worker!(w,3;site_inputs=(:layer,),per_site=(next=Float64,)) do worker
            b=ChunkNL.get_buffers(worker)
            for site in worker
                v=ChunkNL.central_feature(site,Val(:layer))
                for k in 1:ChunkNI.neighbour_count(site)
                    ChunkNI.active(site,k) || continue
                    v+=ChunkNL.neighbor_feature(site,k,Val(:layer))
                end
                b.per_site.next[ChunkNL.site_ordinal(site)]=v
            end
        end
        @test ChunkNL.export_field(w,:next)==[6.,7.,4.,6.]

        ready=Channel{Nothing}(1); release=Channel{Nothing}(1)
        failed_ready=Channel{Nothing}(1)
        failure=ErrorException("deliberate worker failure")
        job=Threads.@spawn try
            ChunkNL.foreach_worker!(w,2;per_site=(failed=Float64,)) do worker
                ChunkNL.get_buffers(worker).per_site.failed[1]=42.
                if first(worker.chunk.centers)==1
                    put!(ready,nothing)
                    take!(release)
                else
                    put!(failed_ready,nothing)
                    throw(failure)
                end
            end
        catch exception
            exception
        end
        take!(ready); take!(failed_ready)
        try
            @test !istaskdone(job)
            @test_throws ArgumentError ChunkNL.foreach_worker!(identity,w,2)
            @test_throws ArgumentError ChunkNL.get_results(w)
            @test_throws ArgumentError close(w)
            @test_throws ArgumentError ChunkNL.prepare_neighbors!(n,sys,1.)
            @test_throws ArgumentError ChunkNL.form_numa_workers!(w,n,3)
            @test_throws ArgumentError ChunkNL.refresh_positions!(w,sys.xyz)
        finally
            put!(release,nothing)
        end
        @test fetch(job)===failure
        @test_throws ArgumentError ChunkNL.get_results(w)
        @test_throws ArgumentError ChunkNL.export_field(w,:failed)
        @test ChunkNL.export_field(w,:layer)==original.layer

        @test_throws ErrorException("deliberate merge failure") ChunkNL.foreach_worker!(
            identity,w,2;scatter=(failed=RejectChunkMerge,),per_site=(rejected=Float64,))
        @test_throws ArgumentError ChunkNL.export_field(w,:rejected)
        @test_throws ArgumentError ChunkNL.foreach_worker!(identity,w,2;per_site=(unassigned=Float64,))
        @test feature_force!(w,2)==original

        X=copy(sys.xyz)
        X[1,1]-=0.01
        old_plan=w.plan
        ChunkNL.refresh_positions!(w,X)
        @test w.plan===old_plan && ChunkNL.export_field(w,:layer)==original.layer
        @test_throws ArgumentError ChunkNL.get_results(w)
        moved=TestSystem(X,sys.H,sys.pbc,sys.labels)
        fresh=ChunkNL.form_numa_workers(ChunkNL.prepare_neighbors(moved,1.;width=4),1;scatter_capacity=3)
        @test feature_force!(w,2)==feature_force!(fresh,1)
        @test feature_force!(w,2).forces!=original.forces
        @test original.layer==[1.,2.,3.,4.]
        n.environment.version+=UInt(1)
        @test_throws ArgumentError ChunkNL.export_field(w,:layer)
        @test_throws ArgumentError ChunkNL.foreach_worker!(identity,w,2)
        @test isempty(w.fields)
    finally
        tasks=isnothing(w.pool) ? Task[] : copy(w.pool.tasks)
        close(w)
        @test all(istaskdone,tasks)
    end

    empty=TestSystem(zeros(3,0),sys.H,sys.pbc,Symbol[])
    no_work=ChunkNL.form_numa_workers(ChunkNL.prepare_neighbors(empty,1.),4)
    ChunkNL.foreach_worker!(x -> error("unexpected work"),no_work,4;
                            sum_outputs=(energy=Float64,),per_site=(x=Float64,))
    @test ChunkNL.get_results(no_work)==(energy=0.,x=Float64[])
    @test isnothing(no_work.pool)
    close(no_work)
end

struct GatedChunkSystem{S} <: AtomsBase.AbstractSystem{3}
    inner::S
    ready::Channel{Nothing}
    release::Channel{Nothing}
    started::Threads.Atomic{Bool}
end
Base.length(s::GatedChunkSystem)=length(s.inner)
AtomsBase.cell_vectors(s::GatedChunkSystem)=AtomsBase.cell_vectors(s.inner)
AtomsBase.periodicity(s::GatedChunkSystem)=AtomsBase.periodicity(s.inner)
AtomsBase.species(s::GatedChunkSystem,i::Integer)=AtomsBase.species(s.inner,i)
function AtomsBase.position(s::GatedChunkSystem,i::Integer)
    if !Threads.atomic_xchg!(s.started,true)
        put!(s.ready,nothing)
        take!(s.release)
    end
    AtomsBase.position(s.inner,i)
end

@testset "Threaded scalar ACE through AtomsCalculators" begin
    model,ps,st=ChunkACE.load_resolved(joinpath(@__DIR__,"../data/ace_ref/resolved_model_v1.toml"))
    sys=TestSystem([6.2 0.4 2.8 7.5;0.3 0.6 1.2 2.5;0.4 1.1 0.2 0.7],
                   18.0Matrix(I,3,3),(true,false,true),[:Al,:Ti,:Al,:Ti])
    tiny=TestSystem([0.2 1.2;0.3 0.65;0.4 0.9],2.0Matrix(I,3,3),(true,true,true),[:Ti,:Al])
    images=TestSystem([0.2 1.2 1.5;0.3 0.65 0.25;0.4 4.7 7.0],diagm([2.0,2.0,20.0]),
                      (true,true,false),[:Ti,:Al,:Ti])
    isolated=TestSystem(reshape([0.2,0.3,0.4],3,1),sys.H,(false,false,false),[:Ti])
    empty=TestSystem(zeros(3,0),sys.H,sys.pbc,Symbol[])
    props=(AC.Energy(),AC.Forces(),AC.Virial())
    function agrees_threaded(y,ref,tol)
        for name in keys(ref)
            name===:state && continue
            @test chunk_numeric(y[name])≈chunk_numeric(ref[name]) atol=tol rtol=tol
        end
    end
    for T in (Float32,Float64), batch_size in (1,3)
        c=ChunkACE.ACECalculator(model,ps,st;precision=T,width=4,slots=2,angular_backend=:block,batch_size)
        state=ChunkACE.ACEEvaluationState(c;scatter_capacity=4)
        tol=T===Float32 ? 3e-4 : 2e-11
        try
            for groups in (1,3), workers in (1,2,4), centers in (nothing,[4,1],Int[])
                r=AC.calculate(props,sys,c,c.parameters,state;executor=:threaded,workers,groups,centers)
                ref=AC.calculate(props,sys,c;executor=:colored,centers)
                @test r.state===state
                @test state.chunks.observed_capacity<=state.scatter_capacity
                agrees_threaded(r,ref,tol)
            end
            for requested in ((AC.Energy(),),(AC.Forces(),),(AC.Virial(),),(AC.Energy(),AC.Virial()))
                r=AC.calculate(requested,sys,c,c.parameters,state;executor=:threaded,workers=3,groups=3)
                ref=AC.calculate(requested,sys,c;executor=:colored)
                agrees_threaded(r,ref,tol)
            end
            for s in (tiny,isolated,empty)
                r=AC.calculate(props,s,c,c.parameters,state;executor=:threaded,workers=3,groups=3)
                agrees_threaded(r,AC.calculate(props,s,c;executor=:colored),tol)
            end
            r=AC.calculate(props,images,c,c.parameters,state;executor=:threaded,workers=3,groups=3,centers=[2,3])
            @test state.chunks.plan.workers>=2
            @test length(state.chunks.groups)>=2
            @test any(i -> length(state.environment.footprints[i]) <
                          count(view(state.environment.active,state.environment.offsets[1,i]:state.environment.offsets[end,i]-1)),
                      (2,3))
            @test !iszero(norm(r.forces[1]))
            agrees_threaded(r,AC.calculate(props,images,c;executor=:colored,centers=[2,3]),tol)
            before=AC.calculate(props,sys,c,c.parameters,state;executor=:threaded,workers=2)
            layout=state.chunks
            pool=layout.pool
            plan=layout.plan
            generation=layout.version
            original=copy(before.forces)
            moved=TestSystem(sys.xyz .+ [0.0 0.01 0.02 0.03;0 0 0 0;0 0 0 0],sys.H,sys.pbc,sys.labels)
            changed=merge(c.parameters,(WB=c.parameters.WB.*T(1.1),))
            r=AC.calculate(props,moved,c,changed,state;executor=:threaded,workers=2,reuse_topology=true)
            @test state.chunks===layout && layout.plan===plan && layout.version==generation
            @test before.forces==original && r.forces!=original
            agrees_threaded(r,AC.calculate(props,moved,c,changed;executor=:colored),tol)
            @test_throws ArgumentError AC.calculate(props,isolated,c,c.parameters,state;
                                                     executor=:threaded,workers=2,reuse_topology=true)
            @test_throws ArgumentError AC.calculate(props,sys,c,c.parameters,state;
                                                     executor=:threaded,workers=2,groups=3,reuse_topology=true)
            # A clean rebuild can change atom mapping and still reuse worker tasks.
            p=reverse(1:length(sys))
            permuted=TestSystem(sys.xyz[:,p],sys.H,sys.pbc,sys.labels[p])
            r=AC.calculate(props,permuted,c,c.parameters,state;executor=:threaded,workers=2)
            @test state.chunks===layout && layout.pool===pool
            @test layout.plan!==plan && layout.version>generation
            @test before.forces==original
            agrees_threaded(r,AC.calculate(props,permuted,c;executor=:colored),tol)
            @test AC.potential_energy(sys,c;executor=:threaded,workers=2)≈before.energy atol=tol*u"eV" rtol=tol
            close(state)
            @test isnothing(state.chunks.pool)
            AC.calculate(props,sys,c,c.parameters,state;executor=:threaded,workers=2)
            @test !isnothing(state.chunks.pool)
        finally
            close(state)
            close(c)
        end
    end

    c=ChunkACE.ACECalculator(model,ps,st;angular_backend=:block,slots=2)
    a=ChunkACE.ACEEvaluationState(c); b=ChunkACE.ACEEvaluationState(c)
    ready=Channel{Nothing}(1); release=Channel{Nothing}(1)
    gated=GatedChunkSystem(sys,ready,release,Threads.Atomic{Bool}(false))
    job=Threads.@spawn AC.calculate(props,gated,c,ps,a;executor=:threaded,workers=2)
    take!(ready)
    try
        @test_throws ArgumentError AC.calculate(props,sys,c,ps,a;executor=:threaded,workers=2)
        @test_throws ArgumentError close(a)
        y=AC.calculate(props,images,c,ps,b;executor=:threaded,workers=2)
        @test b.chunks.plan.workers==2
        agrees_threaded(y,AC.calculate(props,images,c;executor=:colored),2e-11)
    finally
        put!(release,nothing)
        try
            agrees_threaded(fetch(job),AC.calculate(props,sys,c;executor=:colored),2e-11)
        finally
            close(a)
            close(b)
        end
    end
end

@testset "Fresh chunk neighbor and workgroup generations" begin
    # Force a change in the actual cutoff graph with equal atom and chunk counts.
    H=20.0Matrix(I,3,3)
    cell=(true,false,true)
    sys=TestSystem([0.15 0.95 1.75 8.5; 0.2 0.2 0.2 0.2; 0.3 0.3 0.3 0.3],H,cell,[:Ti,:Al,:Ti,:Al])
    moved=TestSystem([19.9 0.5 0.85 8.5; 0.2 0.2 0.2 0.2; 0.3 0.3 0.3 0.3],H,cell,sys.labels)
    small_cell=TestSystem([4.1 0.8; 0.1 0.2; 0.3 0.4],diagm([0.9,0.8,20.]),(true,true,false),[:Al,:Ti])
    n=ChunkNL.prepare_neighbors(sys,1.;width=4)
    w=ChunkNL.form_numa_workers(n,1;scatter_capacity=3)
    run!(layout,workers) = ChunkNL.foreach_worker!(layout,workers;
                scatter=(forces=SVector{3,Float64},),per_site=(degree=Int,)) do worker
        for site in worker
            ChunkNL.get_buffers(worker).per_site.degree[ChunkNL.site_ordinal(site)]=
                count(k -> ChunkNI.active(site,k),1:ChunkNI.neighbour_count(site))
        end
    end
    function images(n,i)
        env=n.environment
        sort!([(env.neighbours[k],Tuple(env.shifts[k])) for k in env.offsets[1,i]:env.offsets[end,i]-1 if env.active[k]])
    end
    initial=images(n,1)
    run!(w,2)
    outputs=ChunkNL.get_results(w)
    original=copy(outputs.degree)
    pool=w.pool
    plan=w.plan
    layout_size=length(w.chunks)
    destinations=first(w.chunks).destinations
    runs=first(w.chunks).runs
    env_storage=(n.environment.neighbours,n.environment.shifts,n.environment.active,n.p,n.inverse)
    try
        @test ChunkNL.prepare_neighbors!(n,moved,1.)===n
        @test images(n,1)!=initial
        @test_throws ArgumentError ChunkNL.export_field(w,:degree)
        @test_throws ArgumentError ChunkNL.foreach_worker!(identity,w,2)
        @test ChunkNL.form_numa_workers!(w,n,1)===w
        @test length(w.chunks)==layout_size
        @test first(w.chunks).destinations===destinations && first(w.chunks).runs===runs
        @test isnothing(w.plan)
        run!(w,2)
        @test w.pool===pool && w.plan!==plan
        @test ChunkNL.get_results(w).degree!=original
        @test maximum(length(buffers[1].forces) for buffers in pool.storage[2])==3
        sparse=TestSystem([0.15 0.60 3.75 8.5; 0.2 0.2 0.2 0.2; 0.3 0.3 0.3 0.3],H,cell,sys.labels)
        ChunkNL.prepare_neighbors!(n,sparse,1.)
        ChunkNL.form_numa_workers!(w,n,1;scatter_capacity=2)
        run!(w,2)
        @test w.pool===pool
        @test maximum(length(buffers[1].forces) for buffers in pool.storage[2])==2
        ChunkNL.form_numa_workers!(w,n,1;scatter_capacity=3)

        invalid=TestSystem(fill(NaN,3,length(sys)),H,cell,sys.labels)
        @test_throws ArgumentError ChunkNL.prepare_neighbors!(n,invalid,1.)
        @test_throws ArgumentError ChunkNL.form_numa_workers!(w,n,1)
        @test_throws ArgumentError ChunkNL.get_results(w)
        for frame in (small_cell,sys,moved)
            @test ChunkNL.prepare_neighbors!(n,frame,1.)===n
            snapshot=ChunkNL.neighbour_list(frame,1.;bruteforce_threshold=32)
            for i in 1:length(frame)
                js,_,ss=ChunkNL.neighbours(snapshot,i)
                @test images(n,i)==sort!([(Int(j),Tuple(Int.(s))) for (j,s) in zip(js,ss)])
            end
            @test !isempty(n.environment.phases)
            ChunkNL.form_numa_workers!(w,n,1)
            run!(w,2)
        end
        @test all(a===b for (a,b) in zip(env_storage,
            (n.environment.neighbours,n.environment.shifts,n.environment.active,n.p,n.inverse)))
        @test_throws ArgumentError ChunkNL.form_numa_workers!(w,n,1;scatter_capacity=1)
        @test_throws ArgumentError ChunkNL.get_results(w)
        ChunkNL.form_numa_workers!(w,n,1;scatter_capacity=3)
        run!(w,2)
        @test outputs.degree==original
    finally
        close(w)
    end
end

@testset "ACE fresh rebuild failure and colored executor switch" begin
    model,ps,st=ChunkACE.load_resolved(joinpath(@__DIR__,"../data/ace_ref/resolved_model_v1.toml"))
    H=20.0Matrix(I,3,3)
    sys=TestSystem([0.25 3.1 8.3; 0.1 0.2 0.3; 0.5 0.7 0.6],H,(true,false,true),[:Ti,:Al,:Ti])
    moved=TestSystem([19.8 3.1 7.9; 0.1 0.2 0.3; 0.5 0.7 0.6],H,sys.pbc,[:Al,:Ti,:Ti])
    c=ChunkACE.ACECalculator(model,ps,st;angular_backend=:block,slots=2)
    state=ChunkACE.ACEEvaluationState(c)
    props=(AC.Energy(),AC.Forces(),AC.Virial())
    try
        previous=AC.calculate(props,sys,c,ps,state;executor=:threaded,workers=2)
        original=copy(previous.forces)
        @test isempty(state.environment.phases)
        state.scatter_capacity=1
        @test_throws ArgumentError AC.calculate(props,moved,c,ps,state;executor=:chunks)
        @test_throws ArgumentError AC.calculate(props,sys,c,ps,state;executor=:threaded,reuse_topology=true)
        @test_throws ArgumentError ChunkNL.get_results(state.chunks)
        state.scatter_capacity=4096
        serial=AC.calculate(props,moved,c,ps,state;executor=:chunks)
        reference=AC.calculate(props,moved,c,ps;executor=:colored,parallel=true)
        @test !isempty(c.environment.phases)
        colored=AC.calculate(props,moved,c,ps,state;executor=:colored,parallel=true)
        @test colored.state===state
        @test_throws ArgumentError AC.calculate(props,moved,c,ps,state;executor=:threaded,reuse_topology=true)
        threaded=AC.calculate(props,moved,c,ps,state;executor=:threaded,workers=2,_inline_merge=true)
        for name in (:energy,:forces,:virial), result in (serial,threaded,colored)
            @test chunk_numeric(result[name])≈chunk_numeric(reference[name]) atol=2e-11 rtol=2e-11
        end
        @test previous.forces==original
    finally
        close(state)
        close(c)
    end
end
