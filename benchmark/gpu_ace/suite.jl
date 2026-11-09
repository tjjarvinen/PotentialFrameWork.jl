# Sequential multi-process orchestration; no GPU/vendor import in the driver.
using TOML
function main()
    root=normpath(joinpath(@__DIR__,"../.."))
    opts=Dict("preset"=>"quick","backend"=>"cpu","data"=>"","output"=>"ace-bench-results","samples"=>"5","seconds"=>"0.3",
        "models"=>"tial,ti,tialv","threads"=>"4","sizes"=>"","audit"=>"true",
        "cpu-project"=>@__DIR__,"amd-project"=>joinpath(@__DIR__,"amd"))
    iseven(length(ARGS)) || error("use --option value pairs")
    for i=1:2:length(ARGS)
        key=replace(ARGS[i],r"^--"=>""); haskey(opts,key) || error("unknown option"); opts[key]=ARGS[i+1]
    end
    opts["preset"] in ("quick","full") || error("preset")
    opts["backend"] in ("cpu","amd","both") || error("backend")
    out=abspath(opts["output"]); mkpath(out); records=Dict[]; omissions=String[]
    jobs=NamedTuple[]
    for backend in (opts["backend"]=="both" ? ("cpu","amd") : (opts["backend"],))
        if backend=="amd" && opts["backend"]=="both"
            project=opts["amd-project"]
            probe=`$(Base.julia_cmd()) --startup-file=no --project=$project -e 'using AMDGPU; exit(AMDGPU.functional() && any(d->occursin("6800 XT",string(d)),AMDGPU.devices()) ? 0 : 2)'`
            if !success(pipeline(addenv(probe,"JULIA_LOAD_PATH"=>"@:@stdlib");stdout=devnull,stderr=devnull))
                push!(omissions,"AMD cases omitted: probe failed in $project; inspect with AMDGPU.versioninfo()")
                continue
            end
        end
        push!(jobs,(;backend,precision="Float64",strategy="staged",threads=opts["threads"],spacing="3.4",sizes=opts["sizes"]))
        if opts["preset"]=="full"
            push!(jobs,(;backend,precision="Float32",strategy="staged",threads=opts["threads"],spacing="3.4",sizes=isempty(opts["sizes"]) ? "4,18" : opts["sizes"]))
            push!(jobs,(;backend,precision="Float64",strategy="staged",threads=opts["threads"],spacing="3.2",sizes=isempty(opts["sizes"]) ? "4,8" : opts["sizes"]))
            if backend=="amd"
                for strategy in ("baseline","compact")
                    push!(jobs,(;backend,precision="Float64",strategy,threads=opts["threads"],spacing="3.4",sizes=isempty(opts["sizes"]) ? "4,8" : opts["sizes"]))
                end
            elseif opts["threads"]!="1"
                push!(jobs,(;backend,precision="Float64",strategy="staged",threads="1",spacing="3.4",sizes=isempty(opts["sizes"]) ? "4,18" : opts["sizes"]))
            end
        end
    end
    for (index,job) in enumerate(jobs)
        id=string(index,"-",job.backend,"-",job.precision,"-",job.strategy,"-t",job.threads,"-a",job.spacing)
        report=joinpath(out,id*".toml"); logfile=joinpath(out,id*".log")
        project=opts[job.backend=="cpu" ? "cpu-project" : "amd-project"]
        args=["--preset",opts["preset"],"--backend",job.backend,"--models",opts["models"],"--data",opts["data"],"--precision",job.precision,
            "--strategy",job.strategy,"--samples",opts["samples"],"--seconds",opts["seconds"],"--output",report,"--spacing",job.spacing]
        isempty(job.sizes) || append!(args,["--sizes",job.sizes])
        cmd=`$(Base.julia_cmd()) --startup-file=no --project=$project --threads=$(job.threads) $(joinpath(@__DIR__,"run.jl")) $args`
        println(id); flush(stdout)
        ok=open(logfile,"w") do io; success(pipeline(addenv(cmd,"JULIA_LOAD_PATH"=>"@:@stdlib");stdout=io,stderr=io)); end
        push!(records,Dict("id"=>id,"passed"=>ok,"report"=>report,"log"=>logfile,"command"=>string(cmd)))
        if ok
            append!(omissions,TOML.parsefile(report)["omitted"])
        end
    end
    if opts["preset"]=="full" && opts["audit"]=="true" && opts["backend"]!="amd"
        for name in split(opts["models"],',')
            path=name=="tial" ? joinpath(root,"data/ace_ref/resolved_model_v1.toml") : joinpath(opts["data"],name*"-resolved.toml")
            isfile(path) || continue
            cmd=`$(Base.julia_cmd()) --startup-file=no --project=$(opts["cpu-project"]) --threads=1 $(joinpath(@__DIR__,"audit.jl")) $path $(joinpath(out,"audit-"*name))`
            log=joinpath(out,"audit-"*name*".log")
            ok=open(log,"w") do io; success(pipeline(addenv(cmd,"JULIA_LOAD_PATH"=>"@:@stdlib");stdout=io,stderr=io)); end
            push!(records,Dict("id"=>"audit-"*name,"passed"=>ok,"log"=>log,"command"=>string(cmd)))
        end
    end
    open(joinpath(out,"suite.toml"),"w") do io
        TOML.print(io,Dict("schema"=>"PotentialFrameWork.GPUACESuite","version"=>1,"options"=>opts,"jobs"=>records,"omitted"=>unique(omissions));sorted=true)
    end
    all(r["passed"] for r in records) || error("one or more benchmark/audit processes failed; inspect suite logs")
end
main()
