using SpheriCart, PotentialFrameWork, BenchmarkTools, InteractiveUtils, TOML, LinearAlgebra, SHA
const A=PotentialFrameWork.ACE
include("workloads.jl")
include("provenance.jl")
include("audit_workload.jl")
function audit(path,out)
    BLAS.set_num_threads(1); mkpath(out)
    model,ps,_=A.load_resolved(path); sys=audit_structure(model.atomic_numbers)
    results=Dict[]
    for T in (Float64,Float32),width in (model.coordinate_kind===:legacy_PI ? (4,) : (1,4,8)),backend in (model.coordinate_kind===:legacy_PI ? (:sphericart,:block) : (:sphericart,:generated,:block)),simd in (false,true),fused in (model.coordinate_kind===:legacy_PI ? (false,) : (false,true))
        width==1 && simd && continue
        fused && !simd && continue
        b=A.adapt_model(model,ps;precision=T)
        env=A.FullEnvironment(;width); input=A.ace_inputs(sys,b.model;environment=env)
        ex=A.CPUExecutor(b.model,env;precision=T,angular_backend=backend)
        X=T.(input.positions); H=T.(input.cell)
        run_case=()->A.execute_prepared!(ex,X,H,b.parameters;simd,fused,parallel=false)
        run_case()
        trial=run(@benchmarkable($run_case(),evals=1);samples=5,seconds=.12)
        push!(results,Dict("precision"=>string(T),"width"=>width,"angular"=>string(backend),"simd"=>simd,"fused"=>fused,
            "atoms"=>length(sys),"median_s"=>BenchmarkTools.median(trial).time*1e-9,"minimum_s"=>minimum(trial).time*1e-9,
            "raw_ns"=>trial.times,"host_bytes"=>trial.memory,"edges"=>count(env.active),"lanes"=>length(env.active),"phase_sizes"=>length.(env.phases)))
        if width==4 && backend===:block && simd && !fused
            s=ex.scratch[1]
            site=A.NL.SiteView(env,X,H,1); siteout=ex.slots[1].outputs[1]; props=(A.AC.Energy(),A.AC.Forces(),A.AC.Virial())
            R=zeros(T,4,length(model.radial_spec)); dR=similar(R)
            jobs=Any[("spline",A.radial_block!,(R,dR,(T(1.1),T(1.7),T(2.3),T(3.1)),b.model.radial[1,1],1,Val(4),(true,true,true,true),A.radial_parameters(b.parameters,:custom_radial))),
                  ("products_reverse",A.product_projection!,(s.P,s.dA,s.dV,s.U,s.V,b.model.products,ex.weights[:,1],b.model.policy.orders,s.polynomial,s.reverse_polynomial)),
                  ("site_math",A.NI.evaluate_site!,(siteout,ex.prepared,site,b.parameters,s,props)),
                  ("force_virial_reduction",A.NL.reduce_site!,(ex.force,ex.inverse,ex.energies,ex.virials,env,1,site,siteout,props)),
                  ("complete_site",A.cpu_site!,(ex,s,1,X,H,b.parameters,Val(width),true,false)),
                  ("density",A.project_rows!,(s,b.model,1,min(4,length(s.Rs)),false)),
                  ("angular",A.angular_block!,(zeros(T,4,(model.lmax+1)^2),zeros(A.SVector{3,T},4,(model.lmax+1)^2),A.BlockAngular(model.lmax;width=4,precision=T),fill(A.SVector{3,T}(1,2,3),4)))]
            # Inspect the actual batched SpheriCart kernels as well as our block path.
            sc=A.angular_basis(model.lmax,T); nr=length(s.Rs)
            norms=zeros(T,nr); unitRs=fill(A.SVector{3,T}(1,2,3)/sqrt(T(14)),nr)
            temps=(x=zeros(T,nr),y=zeros(T,nr),z=zeros(T,nr),r²=zeros(T,nr),
                s=zeros(T,nr,model.lmax+1),c=zeros(T,nr,model.lmax+1),Q=zeros(T,nr,(model.lmax+1)^2),Flm=sc.solids.Flm)
            append!(jobs,[("sphericart_normalize",SpheriCart._normalise_Rs!,(norms,unitRs,sc,s.Rs)),
                ("sphericart_kernel",SpheriCart.solid_harmonics_with_grad!,(s.Y,s.dY,Val(model.lmax),unitRs,temps)),
                ("sphericart_rescale",SpheriCart._rescale_∇Z2∇Y!,(s.dY,unitRs,norms))])
            for (name,f,args) in jobs
                f(args...)
                stem=joinpath(out,name*"-"*string(T))
                open(stem*".llvm.txt","w") do io; code_llvm(io,f,typeof(args);optimize=true); end
                open(stem*".native.txt","w") do io; code_native(io,f,typeof(args)); end
            end
        end
    end
    report=Dict("schema"=>"PotentialFrameWork.CPUSIMDAudit","version"=>1,"julia"=>string(VERSION),"cpu"=>Sys.cpu_info()[1].model,"cpu_target"=>Sys.CPU_NAME,"source_files_sha256"=>source_provenance(normpath(joinpath(@__DIR__,"../.."))),
        "versions"=>Dict(id.name=>string(pkgversion(mod)) for (id,mod) in Base.loaded_modules if id.name in ("SIMD","SpheriCart","BenchmarkTools","StaticArrays","KernelAbstractions","AcceleratedKernels")),
        "threads"=>Threads.nthreads(),"blas_threads"=>BLAS.get_num_threads(),"model_sha256"=>model.source_sha256,
        "commit"=>strip(read(`git rev-parse HEAD`,String)),"cases"=>results,"boundary"=>"prepared numeric EFV; no extraction or output copies; serial executor; warmed")
    open(joinpath(out,"timings.toml"),"w") do io; TOML.print(io,report;sorted=true); end
end
length(ARGS)==2 || error("usage: audit.jl RESOLVED.toml OUTPUT_DIRECTORY")
audit(ARGS...)
