using StaticArrays: SVector
import NeighbourLists as ReferenceNL

function directed_signature(list,i)
    js,Rs,Ss=neighbours(list,i)
    [(Int(js[k]),Tuple(Int.(Ss[k])),Tuple(round.(Float64.(Rs[k]);digits=6))) for k in eachindex(js)] |> Set
end

@testset "directed neighbour interface and image labels" begin
    p=potential()
    cells=([9. 0 0;0 9. 0;0 0 9.],
           [7. 2. 1.;0 6. 1.;0 0 5.],
           [1.7 0.6 0.2;0 1.8 0.4;0 0 1.6],
           [4. 3.8 0.;0 0.6 0;0 0 4.])
    for H in cells, pb in ((true,true,true),(true,false,true),(false,false,false))
        source=fixture(7)
        xyz=H*(source.xyz/9)
        # Noncanonical positions test shifts relative to the ORIGINAL coordinates.
        xyz[:,1]+=3H[:,1]-2H[:,2]
        s=TestSystem(xyz,H,pb,source.labels)
        ref=image_oracle(s,p;shell=12)
        for threshold in (0,32), width in (1,4)
            list=neighbour_list(s,p.radius;width,bruteforce_threshold=threshold)
            raw=neighbour_list([SVector{3}(xyz[:,i]) for i=1:7],p.radius,transpose(H),pb;
                               width,bruteforce_threshold=threshold)
            expected=[Set{Tuple{Int,NTuple{3,Float64}}}() for _=1:7]
            for (i,j,d) in ref.pairs
                push!(expected[i],(j,d)); push!(expected[j],(i,Tuple(-SVector(d))))
            end
            for i=1:7
                js,Rs,Ss=neighbours(list,i)
                @test eltype(js)===Int32
                @test eltype(Rs)===SVector{3,Float64}
                @test eltype(Ss)===SVector{3,Int32}
                actual=Set((Int(js[k]),Tuple(round.(Rs[k];digits=7))) for k in eachindex(js))
                @test actual==expected[i]
                @test length(js)==length(actual) # image pairs must not be duplicated
                @test directed_signature(raw,i)==directed_signature(list,i)
                @test ReferenceNL.neighbours(list,i)==neighbours(list,i)
                for k in eachindex(js)
                    # Keep open axes out of the image sum (including infinite vacuum).
                    reconstructed=xyz[:,js[k]]-xyz[:,i]
                    for a=1:3
                        pb[a] && (reconstructed+=H[:,a]*Ss[k][a])
                        !pb[a] && @test Ss[k][a]==0
                    end
                    @test reconstructed ≈ Rs[k] atol=1e-12
                    @test (Int(i),Tuple(-Ss[k]),Tuple(round.(-Rs[k];digits=6))) in directed_signature(list,js[k])
                end
            end
            @test ReferenceNL.nsites(list)==7
            @test ReferenceNL.npairs(list)==2length(ref.pairs)
            @test_throws BoundsError neighbours(list,0)
            @test_throws BoundsError neighbours(list,8)
        end
    end
end

@testset "reference NeighbourLists convention and snapshot lifetime" begin
    H=[8. 2.2 -0.7;0 7. 1.1;0 0 6.]
    for pb in ((true,true,true),(true,false,true),(false,false,false))
        source=fixture(17)
        xyz=H*(source.xyz/9)
        xyz[:,1]+=3H[:,1]; xyz[:,2]-=2H[:,2]
        X=[SVector{3}(xyz[:,i]) for i=1:17]
        cell=Matrix(transpose(H))
        ours=neighbour_list(X,2.,cell,pb)
        current_api=isdefined(ReferenceNL,:neighbour_list)
        # 0.5 reports shifts relative to its internally wrapped coordinates.
        # Give that legacy reference canonical periodic input, then relabel its
        # shifts with our known windings; current 0.6 is compared directly.
        winding=[SVector(ntuple(a -> !pb[a] ? 0 : i==1 && a==1 ? 3 : i==2 && a==2 ? -2 : 0,3)) for i=1:17]
        reference=if current_api
            ReferenceNL.neighbour_list(X,2.,cell,pb)
        else
            canonical=[X[i]-H*winding[i] for i=1:17]
            ReferenceNL.PairList(canonical,2.,cell,pb)
        end
        for i=1:17
            js,Rs,Ss=current_api ? ReferenceNL.neighbours(reference,i) : ReferenceNL.neigss(reference,i)
            if !current_api
                Ss=[SVector{3,Int}(Ss[k])+winding[i]-winding[js[k]] for k in eachindex(js)]
            end
            expected=Set((Int(js[k]),Tuple(Int.(Ss[k])),Tuple(round.(Rs[k];digits=6))) for k in eachindex(js))
            @test directed_signature(ours,i)==expected
        end
        original=directed_signature(ours,1)
        X[1]+=SVector(9.,-4.,2.); cell[1,1]+=1
        @test directed_signature(ours,1)==original
        js,Rs,Ss=neighbours(ours,1)
        empty!(js); empty!(Rs); empty!(Ss)
        @test directed_signature(ours,1)==original
    end
    # Matrix positions, explicit units and Float32/index-type selection.
    X=Float32[0 0.5;0 0.2;0 0.1]
    cell=Float32[5 0 0;0 5 0;0 0 5]
    list=neighbour_list(X*u"nm",1u"nm",cell*u"nm",(false,false,false);precision=Float32,int_type=Int64)
    js,Rs,Ss=neighbours(list,1)
    @test js==[2]
    @test eltype(Rs)===SVector{3,Float32}
    @test eltype(Ss)===SVector{3,Int64}
    @test Rs[1] ≈ X[:,2]*10
    # The same atom appears as both directed images; no zero-image self entry.
    s=TestSystem(zeros(3,1),[1.2 0 0;0 0 0;0 0 Inf],(true,false,false),[:H])
    list=neighbour_list(s,2.)
    js,Rs,Ss=neighbours(list,1)
    @test js==[1,1]
    @test Set(Tuple.(Ss))==Set(((1,0,0),(-1,0,0)))
    @test sort(first.(Rs)) ≈ [-1.2,1.2]
    @test all(r -> iszero(r[2]) && iszero(r[3]),Rs)
end

@testset "neighbour interface validation" begin
    cell=[8. 0 0;0 8. 0;0 0 8.]; pb=(true,true,true)
    empty=neighbour_list(SVector{3,Float64}[],2.,cell,pb)
    @test ReferenceNL.nsites(empty)==0
    @test ReferenceNL.npairs(empty)==0
    @test_throws BoundsError neighbours(empty,1)
    @test isempty(first(neighbours(neighbour_list(zeros(3,1),2.,cell,pb),1)))
    @test isempty(first(neighbours(neighbour_list([SVector(0.,0.,0.),SVector(2.,0.,0.)],2.,cell,pb),1)))
    @test_throws DimensionMismatch neighbour_list(zeros(2,3),2.,cell,pb)
    @test_throws DimensionMismatch neighbour_list([[0.,0.]],2.,cell,pb)
    @test_throws DimensionMismatch neighbour_list(zeros(3,0),2.,zeros(2,2),pb)
    @test_throws DimensionMismatch neighbour_list(zeros(3,0),2.,cell,(true,true))
    @test_throws ArgumentError neighbour_list(zeros(3,0),2.,cell,pb;int_type=UInt32)
    @test_throws ArgumentError neighbour_list(zeros(3,128),2.,cell,pb;int_type=Int8)
    @test_throws ArgumentError neighbour_list(zeros(3,0),1u"eV",cell,pb)
    @test_throws ArgumentError neighbour_list(zeros(3,0),0.,cell,pb)
    # Atom indices fit Int8, but an unwrapped position requires an overflowing shift.
    @test_throws InexactError neighbour_list([SVector(0.,0.,0.),SVector(1600.1,0.,0.)],2.,cell,pb;int_type=Int8)
    # -128 fits Int8 but its reverse +128 does not; negation must not wrap.
    @test_throws OverflowError neighbour_list([SVector(0.,0.,0.),SVector(1024.1,0.,0.)],2.,cell,pb;int_type=Int8)
end
