@testset "user neighbour-strategy benchmark" begin
    H=[1.7 0.6 0.2;0 1.8 0.4;0 0 1.6]
    source=fixture(5)
    skew=TestSystem(H*(source.xyz/9),H,source.pbc,source.labels)
    systems=(fixture(0),fixture(1),fixture(7),fixture(40;pbc=(false,false,false)),skew)
    for s in systems, parallel in (false,true)
        original=copy(s.xyz); cell=copy(s.H); labels=copy(s.labels)
        io=IOBuffer()
        r=benchmark_neighbourlist(s,0.2u"nm";samples=2,warmup=1,parallel,io)
        @test r.pairs==length(image_oracle(s,potential();shell=3).pairs)
        @test r.metadata.atoms==length(s)
        @test r.metadata.cutoff_angstrom==2.
        @test r.metadata.parallel==parallel
        @test r.metadata.threads==Threads.nthreads()
        @test r.metadata.samples==2
        @test r.preferred in (:bruteforce,:cells,:equivalent)
        @test r.cells_over_bruteforce==r.cells.median/r.bruteforce.median
        for method in (r.bruteforce,r.cells)
            @test length(method.seconds)==length(method.bytes)==2
            @test 0<method.minimum<=method.median<=method.maximum
            @test method.preparation_median>=0
            @test method.traversal_median>=0
            @test method.bytes_median>=0
            @test all(>=(0),method.bytes)
        end
        output=String(take!(io))
        @test occursin("bruteforce",output) && occursin("cells",output)
        @test occursin("Excludes potential evaluation",output)
        @test s.xyz==original && s.H==cell && isequal(s.labels,labels)
    end
    r=benchmark_neighbourlist(skew,2.;width=1,precision=Float32,samples=1,warmup=1,io=nothing)
    @test r.metadata.width==1
    @test r.metadata.precision===Float32
    @test r.pairs==length(image_oracle(skew,potential();shell=3).pairs)
    s=fixture(3)
    for kw in ((;samples=0),(;samples=1.5),(;warmup=0),(;warmup=1.5),(;width=3))
        @test_throws ArgumentError benchmark_neighbourlist(s,2.;kw...,io=nothing)
    end
    for cutoff in (0.,-1.,Inf,NaN,1u"eV")
        @test_throws ArgumentError benchmark_neighbourlist(s,cutoff;io=nothing)
    end
end
