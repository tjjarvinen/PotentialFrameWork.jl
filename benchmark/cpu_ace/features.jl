include("common.jl")
using ChainRulesCore, ForwardDiff
function main()
    BLAS.set_num_threads(1)
    output=ARGS[1]; samples=length(ARGS)>1 ? parse(Int,ARGS[2]) : 10
    seconds=length(ARGS)>2 ? parse(Float64,ARGS[3]) : 3.
    radial=A.CallableRadial(A.SmoothExponential(),
        A.RadialMetadata([(1,0),(1,1),(2,1)];cutoff=4.,ad=true);
        config=(cutoff=4.,powers=[0,0,1]))
    layer=A.ACEFeatureLayer(species=[:Si,:O],radial=radial,correlation=2,lmax=1,
        output_irreps=((L=0,parity=:even,channels=3),(L=1,parity=:odd,channels=2),
                      (L=1,parity=:even,channels=2),(L=2,parity=:even,channels=2)))
    ps=A.feature_parameters(layer;custom_radial=(rate=[.4],amplitude=[1.3]))
    records=Dict[]
    for n in (3,32)
        X=[.1+.67*((i-1)%4)+.03sin(i*a) for a=1:3,i=1:n]
        # Noncollinear bounded cluster, functional feature path only.
        X[2,:].=[.2+.63*((i-1)÷4%4) for i=1:n]
        X[3,:].=[.3+.71*((i-1)÷16) for i=1:n]
        sys=AB.isolated_system([[:Si,:O][mod1(i,2)]=>A.SVector{3}(X[:,i])*U.u"Å" for i=1:n])
        input=A.ace_inputs(sys,layer)
        f=()->A.ace_features(input.positions,layer,ps,input.environment,input.cell)
        rule=()->rrule(A.ace_features,input.positions,layer,ps,input.environment,input.cell)
        y,pb=rule(); cot=map(b->reshape(sin.(1:length(b)),size(b)),y)
        @assert all(b->all(isfinite,b),y)
        push!(records,Dict("atoms"=>n,"density_count"=>length(layer.density),"product_count"=>length(layer.products),
            "irreps"=>string(layer.irreps),"forward"=>measure(f,samples,seconds),
            "rule_forward"=>measure(rule,samples,seconds),"position_parameter_pullback"=>measure(()->pb(cot),samples,seconds),
            "output_payload_bytes"=>sum(sizeof,y),"saved_pullback_summarysize"=>Base.summarysize(pb),
            "boundary"=>"functional serial CPU features, owning outputs; prepared environment; joint position/parameter pullback"))
    end
    open(io->TOML.print(io,Dict("schema"=>"PFW.CPUACE.Features.v1","julia"=>string(VERSION),"threads"=>Threads.nthreads(),"blas_threads"=>BLAS.get_num_threads(),"cases"=>records);sorted=true),output,"w")
end
main()
