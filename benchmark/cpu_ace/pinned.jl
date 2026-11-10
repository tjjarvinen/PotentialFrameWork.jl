# Targeted controlled-affinity repeat for observed parallel baseline variation.
# taskset -c 0-3 julia --project=ENV -t4 pinned.jl current MODEL_CACHE ORACLE_CACHE OUTPUT
include("common.jl")
function main()
    mode,data,oracles,output=ARGS
    mode in ("current","historical") || error("mode")
    BLAS.set_num_threads(1); records=Dict[]
    for name in ("small-tial","medium-ti","large-tial")
        path=joinpath(data,name*"-resolved.toml"); m,ps,_=A.load_resolved(path)
        d=bcc_geometry(8,m.atomic_numbers); sys=structure(d)
        env=A.FullEnvironment(width=4); input=A.ace_inputs(sys,m;environment=env)
        w=A.CPUExecutor(m,env)
        oracle=TOML.parsefile(joinpath(oracles,name*"-oracle.toml"))
        raw=only(filter(r->r["g"]==8,oracle["records"]))
        ref=(;energy=raw["energy"],forces=reshape(raw["forces"],3,:),virial=reshape(raw["virial"],3,3))
        for algorithm in (mode=="historical" ? (:historical,) : (:reference,:optimized))
            kw=mode=="historical" ? (;parallel=true) : (;parallel=true,product_algorithm=algorithm)
            f=()->A.execute_prepared!(w,input.positions,input.cell,ps;kw...)
            f(); y=A.prepared_result(w); agree(y,ref,Float64)
            calc=mode=="historical" ? A.ACECalculator(m,ps) : A.ACECalculator(m,ps;product_algorithm=algorithm)
            push!(records,Dict("model"=>name,"atoms"=>1024,"precision"=>"Float64","algorithm"=>string(algorithm),
                "prepared"=>measure(f,10,60.),"complete_public"=>measure(()->AC.energy_forces_virial(sys,calc;parallel=true),10,60.),
                "errors"=>accuracy(y,ref),"metadata"=>metadata(m,ps),"model_sha256"=>bytes2hex(sha256(read(path)))))
        end
    end
    report=Dict("schema"=>"PFW.CPUACE.PinnedRepeat.v1","threads"=>Threads.nthreads(),"blas_threads"=>BLAS.get_num_threads(),
        "affinity"=>strip(read(`taskset -pc $(getpid())`,String)),"source_files_sha256"=>source_provenance(dirname(dirname(pathof(PotentialFrameWork)))),
        "script_sha256"=>bytes2hex(sha256(read(@__FILE__))),"common_sha256"=>bytes2hex(sha256(read(joinpath(@__DIR__,"common.jl")))),
        "loaded_package_source"=>pathof(PotentialFrameWork),"cases"=>records,"julia"=>string(VERSION),"cpu"=>Sys.cpu_info()[1].model)
    mkpath(dirname(output)); open(io->TOML.print(io,report;sorted=true),output,"w")
end
main()
