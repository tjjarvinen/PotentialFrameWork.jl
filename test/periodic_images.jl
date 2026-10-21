import AtomsCalculators

# Independent exhaustive image sum: fixed, deliberately oversized lattice cube;
# original atom order, no reciprocal bounds, sorting, cell stencil or SIMD masks.
function image_oracle(s,p; shell=10)
    n=length(s); f=zeros(3,n); v=zeros(3,3); energy=0.0
    pairs=Set{Tuple{Int,Int,NTuple{3,Float64}}}()
    axes=ntuple(a -> s.pbc[a] ? (-shell:shell) : (0:0),3)
    for i=1:n, j=i:n, z in axes[3], y in axes[2], x in axes[1]
        shift=(x,y,z)
        i==j && shift<=(0,0,0) && continue
        # Skip open vectors entirely, so zero/infinite vacuum vectors are valid.
        d=s.xyz[:,j]-s.xyz[:,i]
        for a=1:3
            s.pbc[a] && (d+=shift[a]*s.H[:,a])
        end
        r2=sum(abs2,d)
        r2<p.radius^2 || continue
        push!(pairs,(i,j,Tuple(round.(d;digits=7))))
        a=p.amplitudes[findfirst(isequal(s.labels[i]),p.labels),findfirst(isequal(s.labels[j]),p.labels)]
        q=1-r2/p.radius^2
        energy+=a*q^3
        force=(-6*a*q^2/p.radius^2)*d
        f[:,i]+=force; f[:,j]-=force
        v-=d*transpose(force)
    end
    (;energy,forces=f,virial=v,pairs)
end

function image_ownership(w::Workspace{W}) where W
    counts=Dict{Tuple{Int,Int,NTuple{3,Float64}},Int}()
    saved=copy(w.force)
    bound=bind_species(potential(),w.species_map)
    for task in eachindex(w.ranges)
        actual=Set([task])
        for i in w.ranges[task], image in w.images[task]
            r=w.ranges[image.cell]
            for j=max(first(r),i+(image.self ? 0 : 1)):W:last(r)
                js,b=NL.geometry(w,i,j,last(r),image)
                for k=1:W
                    NL.lane(b.active,k) || continue
                    a,bb=w.p[i],w.p[js[k]]
                    d=[NL.lane(b.dx,k),NL.lane(b.dy,k),NL.lane(b.dz,k)]
                    if a>bb; a,bb=bb,a; d=-d; end
                    key=(a,bb,Tuple(round.(d;digits=7)))
                    counts[key]=get(counts,key,0)+1
                    push!(actual,image.cell)
                end
            end
        end
        @test issubset(actual,Set(w.footprints[task]))
        fill!(w.force,0)
        PP.cell_task!(w,bound,EnergyForcesVirial(),task)
        destinations=Set(c for c in eachindex(w.ranges) if any(!iszero,w.force[:,w.ranges[c]]))
        @test issubset(destinations,Set(w.footprints[task]))
    end
    w.force.=saved
    counts
end

@testset "general cells, explicit images and brute-force dispatch" begin
    p=potential()
    cells=(
        [7. 2.5 -1.; 0. 6. 1.5; 0. 0. 5.], # ordinary triclinic
        [4. 3.8 0.; 0. 0.6 0.; 0. 0. 4.], # rounding fractional components is insufficient
        [1.2 0.4 0.; 0. 1.3 0.3; 0. 0. 1.1], # multiple images incl. self images
        [2. 0 0; 0 2. 0; 0 0 2.], # exact self-image cutoff
        [0. -5. 0.; 6. 0. 1.; 0. 0. 7.], # rotated/left handed generators
    )
    for H in cells, pbc in ((true,true,true),(true,false,true),(false,true,false))
        frac=[mod(sin(i*(a+0.37))*31.3,1) for a=1:3,i=1:7]
        s=TestSystem(H*frac,H,pbc,[:H,:He,:H,:He,:H,:He,:H])
        ref=image_oracle(s,p)
        @test ref.pairs==image_oracle(s,p;shell=12).pairs
        for threshold in (0,32), width in (1,4,8)
            w=Workspace(;width,bruteforce_threshold=threshold)
            for parallel in (false,true)
                agrees(evaluate!(w,s,p;parallel),ref;tol=2e-10)
                agrees(evaluate!(w,s,p;parallel,request=EnergyOnly()),ref;tol=2e-10)
                @test validate_schedule(w)
            end
            @test w.algorithm==(threshold==0 ? :cells : :bruteforce)
            @test sizeof(w.force)==3*length(s)*sizeof(Float64)
            counts=image_ownership(w)
            @test Set(keys(counts))==ref.pairs
            @test all(==(1),values(counts))
            ps=AtomsCalculators.get_parameters(p)
            @test parameter_energy(s,p,ps;neighbourlist=w,parallel=true) ≈ ref.energy atol=2e-10
        end
        translated=TestSystem(s.xyz .+ H[:,findfirst(pbc)]*2,H,pbc,s.labels)
        agrees(evaluate(translated,p),ref;tol=2e-10)
        perm=reverse(1:length(s))
        rp=evaluate(TestSystem(s.xyz[:,perm],H,pbc,s.labels),p)
        @test rp.energy ≈ ref.energy*u"eV" atol=2e-10*u"eV"
        @test ustrip.(u"eV/Å",rp.forces) ≈ ref.forces[:,perm] atol=2e-10
    end
    # Explicit one-atom periodic chain: ±a are one self-image pair, zero net force.
    H=[1.2 0 0;0 0 0;0 0 0.]
    s=TestSystem(zeros(3,1),H,(true,false,false),[:H])
    r=evaluate(s,p)
    @test r.energy ≈ (1-(1.2/2)^2)^3*u"eV"
    @test all(iszero,r.forces)
    @test r.virial[1,1]>0u"eV"
    agrees(r,image_oracle(s,p))
    # Slab with no meaningful open-axis cell vector.
    H=[3. 1. 0;0 2.5 0;0 0 Inf]
    s=TestSystem([0. 0.5;0 0.2;0 0.8],H,(true,true,false),[:H,:He])
    agrees(evaluate(s,p),image_oracle(s,p))
    # Nonperiodic and minimum-image orthorhombic brute/cell paths agree.
    for pbc in ((false,false,false),(true,false,true),(true,true,true))
        s=fixture(17;pbc)
        for threshold in (0,16,17,32), width in (1,4)
            w=Workspace(;width,bruteforce_threshold=threshold)
            agrees(evaluate!(w,s,p;parallel=true),oracle(s,p))
            @test w.algorithm==(length(s)<=threshold ? :bruteforce : :cells)
        end
    end
    for T in (Float32,Float64), n in (0,1,3,35), threshold in (0,32)
        H=[1.2 0.4 0.;0 1.3 0.3;0 0 1.1]
        s=TestSystem(H*[mod(sin(i*a)*17.3,1) for a=1:3,i=1:n],H,(true,true,true),fill(:H,n))
        agrees(evaluate(s,p;precision=T,bruteforce_threshold=threshold),image_oracle(s,p);tol=T==Float32 ? 2e-2 : 1e-9)
    end
    @test_throws ArgumentError Workspace(bruteforce_threshold=-1)
end

@testset "skew multi-image force and strain derivatives" begin
    H=[1.7 0.6 0.2;0 1.8 0.4;0 0 1.6]
    s=TestSystem(H*[0.1 0.4 0.8;0.2 0.6 0.3;0.1 0.7 0.5],H,(true,true,true),[:H,:He,:H])
    p=potential(); h=1e-5
    for threshold in (0,32)
        r=evaluate(s,p;bruteforce_threshold=threshold)
        for i=1:length(s), a=1:3
            xp=copy(s.xyz); xm=copy(s.xyz); xp[a,i]+=h; xm[a,i]-=h
            ep=evaluate(TestSystem(xp,H,s.pbc,s.labels),p;request=EnergyOnly(),bruteforce_threshold=threshold).energy
            em=evaluate(TestSystem(xm,H,s.pbc,s.labels),p;request=EnergyOnly(),bruteforce_threshold=threshold).energy
            @test ustrip(r.forces[a,i]) ≈ -ustrip(ep-em)/(2h) atol=2e-8
        end
        I3=[1. 0 0;0 1. 0;0 0 1.]
        for strain in (I3,[0. 1. 0;0 0 0;0 0 0.],[0. 0 0;0 0 1.;0 0 0])
            ap=I3+h*strain; am=I3-h*strain
            ep=evaluate(TestSystem(ap*s.xyz,ap*H,s.pbc,s.labels),p;request=EnergyOnly(),bruteforce_threshold=threshold).energy
            em=evaluate(TestSystem(am*s.xyz,am*H,s.pbc,s.labels),p;request=EnergyOnly(),bruteforce_threshold=threshold).energy
            @test ustrip(ep-em)/(2h) ≈ -sum(ustrip.(r.virial).*transpose(strain)) atol=3e-8
        end
    end
end

@testset "skew-cell scheduling concurrency and workspace path changes" begin
    s0=fixture(257;L=30.,cluster=true)
    H=[30. 8. 3.;0 28. 6.;0 0 25.]
    s=TestSystem(H*(s0.xyz/30),H,(true,true,true),s0.labels)
    p=potential(); w=Workspace(bruteforce_threshold=0)
    ref=image_oracle(s,p;shell=1)
    for repeat=1:4
        agrees(evaluate!(w,s,p;parallel=true),ref)
        @test validate_schedule(w)
        @test maximum(length,w.phases)>1
    end
    # Reusing a workspace must clear both minimum-image and explicit-image metadata.
    for target in (fixture(3),s,fixture(1;L=1.2),fixture(0),fixture(45))
        agrees(evaluate!(w,target,p;parallel=true),image_oracle(target,p;shell=2);tol=2e-10)
        @test validate_schedule(w)
    end
    s=fixture(7)
    w=Workspace(bruteforce_threshold=0)
    @test AtomsCalculators.potential_energy(s,p;bruteforce_threshold=0) ≈ evaluate(s,p).energy
    @test AtomsCalculators.potential_energy(s,p;neighbourlist=w,bruteforce_threshold=0) ≈ evaluate(s,p).energy
    @test_throws ArgumentError AtomsCalculators.potential_energy(s,p;neighbourlist=w,bruteforce_threshold=1)
    @test_throws ArgumentError parameter_energy(s,p,AtomsCalculators.get_parameters(p);neighbourlist=w,bruteforce_threshold=1)
end
