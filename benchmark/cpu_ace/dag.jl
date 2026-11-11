# CPU DAG controlled comparison. Run jobs sequentially, with taskset externally.
include("common.jl")
using Profile
function dag_metadata(m,plan,T,w)
    n=plan.nleaves; nodes=length(plan.left); depth=zeros(Int,n+nodes)
    for i=1:nodes; depth[n+i]=1+max(depth[plan.left[i]],depth[plan.right[i]]); end
    terminals=Set(plan.terminals)
    factors=sum(length,m.products;init=0)
    payload=sum(sizeof,(plan.left,plan.right,plan.terminals))+sum(sizeof,plan.topology)
    scratch=w.scratch[1]
    Dict("leaves"=>n,"unique_terminals"=>length(terminals),"multiplication_nodes"=>nodes,
        "intermediate_only_nodes"=>count(i->!(n+i in terminals),1:nodes),"max_depth"=>maximum(depth;init=0),
        "original_factor_count"=>factors,"prefix_forward_multiplications"=>factors,
        "prefix_reverse_multiplications"=>2factors,"dag_forward_multiplications"=>nodes,
        "dag_reverse_multiplications"=>2nodes,"prefix_reverse_additions"=>2factors,
        "dag_reverse_additions"=>2nodes,"dag_seed_additions"=>length(m.products),
        "dag_leaf_additions"=>n,"prefix_leaf_additions"=>2factors,"dag_clear_elements"=>n+nodes,"prefix_clear_elements"=>2n,
        "sharing_ratio_vs_independent_nontrivial_products"=>(sum(p->max(length(p)-1,0),m.products;init=0)/max(nodes,1)),
        "plan_array_payload_bytes"=>payload,"dag_per_slot_bytes"=>sizeof(scratch.dag_values)+sizeof(scratch.dag_adjoints),
        "retained_polynomial_per_slot_bytes"=>sizeof(scratch.polynomial)+sizeof(scratch.reverse_polynomial),
        "total_per_slot_bytes"=>sum(sizeof(getfield(scratch,k)) for k in fieldnames(typeof(scratch))),
        "all_slot_bytes"=>sum(s->sum(sizeof(getfield(s,k)) for k in fieldnames(typeof(s))),w.scratch),
        "single_shared_force_bytes"=>sizeof(w.force))
end
function main()
    opts=Dict("models"=>"small-tial,medium-ti,large-tial,ti,tialv","sizes"=>"4,8","precisions"=>"Float64,Float32",
        "algorithms"=>"reference,optimized,dag","samples"=>"10","seconds"=>"120","diagnostics"=>"true",
        "profile"=>"false","codegen"=>"false","output"=>"dag.toml","stageone"=>"false","parallel"=>"false")
    iseven(length(ARGS)) || error("--key value pairs")
    for i=1:2:length(ARGS); k=replace(ARGS[i],"--"=>""); haskey(opts,k)||error(k); opts[k]=ARGS[i+1]; end
    BLAS.set_num_threads(1); samples=parse(Int,opts["samples"]); seconds=parse(Float64,opts["seconds"])
    out=abspath(opts["output"]); mkpath(dirname(out)); records=Dict[]
    root=dirname(dirname(pathof(PotentialFrameWork)))
    report=Dict("schema"=>"PFW.CPUACE.DAG.v1","options"=>opts,"julia"=>string(VERSION),
        "cpu"=>Sys.cpu_info()[1].model,"cpu_target"=>Sys.CPU_NAME,"llvm"=>string(Base.libllvm_version),
        "threads"=>Threads.nthreads(),"blas_threads"=>BLAS.get_num_threads(),"affinity"=>strip(read(`taskset -pc $(getpid())`,String)),
        "source_files_sha256"=>source_provenance(root),"harness_sha256"=>Dict(f=>bytes2hex(sha256(read(joinpath(@__DIR__,f)))) for f in ("dag.jl","common.jl")),
        "versions"=>Dict(id.name=>string(pkgversion(mod)) for (id,mod) in Base.loaded_modules if id.name in ("SIMD","SpheriCart","BenchmarkTools","StaticArrays","KernelAbstractions","AcceleratedKernels","AtomsBase","AtomsCalculators")),"cases"=>records)
    save()=open(io->TOML.print(io,report;sorted=true),out,"w")
    for name in split(opts["models"],','), g in parse.(Int,split(opts["sizes"],',')), tname in split(opts["precisions"],',')
        T=tname=="Float32" ? Float32 : Float64
        path=joinpath(@__DIR__,"results/2026-11-10/models",name*"-resolved.toml")
        m,ps,_=A.load_resolved(path;precision=T)
        d=bcc_geometry(g,m.atomic_numbers); sys=structure(d;T)
        env=A.FullEnvironment(width=4); input=A.ace_inputs(sys,m;environment=env,precision=T)
        X=T.(input.positions); H=T.(input.cell)
        w=A.CPUExecutor(m,env;precision=T)
        A.execute_prepared!(w,X,H,ps;product_algorithm=:reference); reference=A.prepared_result(w)
        oraclepath=joinpath(@__DIR__,"results/2026-11-10/oracles-final",name*"-oracle.toml")
        oracle=TOML.parsefile(oraclepath); oi=findfirst(r->r["g"]==g,oracle["records"])
        external=isnothing(oi) ? nothing : let r=oracle["records"][oi]
            (;energy=r["energy"],forces=reshape(r["forces"],3,:),virial=reshape(r["virial"],3,3))
        end
        for algorithm in Symbol.(split(opts["algorithms"],','))
            f=()->A.execute_prepared!(w,X,H,ps;product_algorithm=algorithm,parallel=opts["parallel"]=="true")
            f(); y=A.prepared_result(w); agree(y,reference,T)
            record=Dict("model"=>name,"atoms"=>length(sys),"g"=>g,"precision"=>tname,"width"=>4,
                "algorithm"=>string(algorithm),"effective_algorithm"=>string(algorithm),"simd"=>true,"fusion"=>false,"angular"=>"sphericart",
                "model_sha256"=>bytes2hex(sha256(read(path))),"geometry_sha256"=>bytes2hex(sha256(reinterpret(UInt8,vec(X)))),
                "cell_sha256"=>bytes2hex(sha256(reinterpret(UInt8,vec(H)))),"species_sha256"=>bytes2hex(sha256(string(d.species))),
                "metadata"=>metadata(m,ps),"reference_errors"=>accuracy(y,reference),"prepared"=>measure(f,samples,seconds))
            if !isnothing(external)
                agree(y,external,T); record["external_errors"]=accuracy(y,external)
                raw=oracle["records"][oi]; record["external_cached_raw_s"]=raw["raw_ns"].*1e-9
                record["external_version"]=oracle["ACEpotentials"]; record["external_boundary"]=raw["boundary"]
                record["oracle_file_sha256"]=bytes2hex(sha256(read(oraclepath)))
            end
            calc=A.ACECalculator(m,ps;precision=T,width=4,product_algorithm=algorithm)
            public=()->AC.energy_forces_virial(sys,calc;parallel=opts["parallel"]=="true")
            py=public(); yn=(;energy=U.ustrip(U.u"eV",py.energy),forces=reduce(hcat,[U.ustrip.(U.u"eV/Å",v) for v in py.forces]),virial=Matrix(U.ustrip.(U.u"eV",py.virial)))
            agree(yn,y,T); record["public_prepared_errors"]=accuracy(yn,y)
            record["complete_public"]=measure(public,samples,seconds)
            fresh=()->AC.energy_forces_virial(sys,A.ACECalculator(m,ps;precision=T,width=4,product_algorithm=algorithm);parallel=opts["parallel"]=="true")
            record["construct_and_public"]=measure(fresh,samples,seconds)
            if opts["stageone"]!="true"
                record["effective_algorithm"]=string(A.effective_product_algorithm(algorithm,m.policy.orders))
                record["graph"]=dag_metadata(m,w.prepared.product_plan,T,w)
                if g==first(parse.(Int,split(opts["sizes"],','))) && algorithm===:dag
                    record["plan_build"]=measure(()->A.ProductDAG(m.products,length(m.density)),samples,seconds)
                    record["executor_build"]=measure(()->A.CPUExecutor(m,env;precision=T,product_algorithm=:dag),samples,seconds)
                end
                if opts["diagnostics"]=="true" && g==4
                    codepath=opts["codegen"]=="true" && algorithm===:dag ? joinpath(dirname(out),"codegen",name*"-"*tname*"-"*string(algorithm)) : ""
                    record["stages"]=diagnostics(w,X,H,ps,algorithm,samples,seconds,codepath)
                end
                if opts["profile"]=="true" && algorithm in (:optimized,:dag)
                    for (boundary,fn) in (("prepared",f),("public",public))
                        Profile.clear(); Profile.@profile for _=1:30; fn(); end
                        open(joinpath(dirname(out),name*"-"*tname*"-"*string(length(sys))*"-"*string(algorithm)*"-"*boundary*"-t"*string(Threads.nthreads())*".profile.txt"),"w") do io
                            Profile.print(io;format=:flat,sortedby=:count)
                        end
                    end
                end
            end
            push!(records,record); save()
            println(name," ",length(sys)," ",tname," ",algorithm," prepared=",record["prepared"]["median_s"]," public=",record["complete_public"]["median_s"]); flush(stdout)
        end
    end
end
main()
