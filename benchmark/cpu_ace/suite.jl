# Sequential process driver, including representative thread scaling.
using TOML
function main()
    opts=Dict("preset"=>"quick","project"=>@__DIR__,"data"=>joinpath(@__DIR__,"results/2026-11-10/models"),
        "output"=>"cpu-ace-suite","oracle-data"=>"","samples"=>"10","seconds"=>"3","sizes"=>"","models"=>"small-tial,medium-ti,large-tial",
        "threads"=>"","widths"=>"4","precisions"=>"Float64,Float32","profile"=>"false","baseline-project"=>"","codegen"=>"false","features"=>"true")
    iseven(length(ARGS)) || error("use --option value pairs")
    for i=1:2:length(ARGS)
        k=replace(ARGS[i],r"^--"=>""); haskey(opts,k) || error("unknown option $k"); opts[k]=ARGS[i+1]
    end
    opts["preset"] in ("quick","full") || error("preset")
    out=abspath(opts["output"]); mkpath(out); jobs=Dict[]
    ts=isempty(opts["threads"]) ? (opts["preset"]=="quick" ? [1] : [1,4,16]) : parse.(Int,split(opts["threads"],','))
    all(>(0),ts) || error("thread counts must be positive")
    # Explicit serial first; no concurrent expensive cases.
    ts=sort(unique(vcat(1,ts)))
    for (label,project,historical) in (("current",opts["project"],false),("historical",opts["baseline-project"],true))
        isempty(project) && continue
        for t in ts
            id=label*"-t"*string(t); file=joinpath(out,id*".toml"); log=joinpath(out,id*".log")
            arguments=["--preset",opts["preset"],"--data",abspath(opts["data"]),"--models",opts["models"],"--samples",opts["samples"],"--seconds",opts["seconds"],
                "--output",file,"--parallel",string(t>1),"--historical",string(historical),"--diagnostics",string(t==1 && !historical),"--widths",opts["widths"],"--precisions",opts["precisions"]]
            isempty(opts["oracle-data"]) || append!(arguments,["--oracle-data",abspath(opts["oracle-data"])])
            isempty(opts["sizes"]) || append!(arguments,["--sizes",opts["sizes"]])
            opts["codegen"]=="true" && !historical && t==1 && append!(arguments,["--codegen",joinpath(out,"codegen")])
            opts["profile"]=="true" && !historical && t==1 && append!(arguments,["--profile",joinpath(out,"profiles")])
            cmd=`$(Base.julia_cmd()) --startup-file=no --project=$project --threads=$t $(joinpath(@__DIR__,"run.jl")) $arguments`
            println(id); flush(stdout)
            started=time()
            passed=open(log,"w") do io; success(pipeline(cmd;stdout=io,stderr=io)); end
            push!(jobs,Dict("id"=>id,"passed"=>passed,"command"=>string(cmd),"report"=>file,"log"=>log,"process_wall_s"=>time()-started))
            open(io->TOML.print(io,Dict("options"=>opts,"jobs"=>jobs);sorted=true),joinpath(out,"suite.toml"),"w")
        end
    end
    if opts["features"]=="true"
        cmd=`$(Base.julia_cmd()) --startup-file=no --project=$(opts["project"]) --threads=1 $(joinpath(@__DIR__,"features.jl")) $(joinpath(out,"features.toml")) $(opts["samples"]) $(opts["seconds"])`
        passed=open(joinpath(out,"features.log"),"w") do io; success(pipeline(cmd;stdout=io,stderr=io)); end
        push!(jobs,Dict("id"=>"features","passed"=>passed,"command"=>string(cmd)))
    end
    open(io->TOML.print(io,Dict("options"=>opts,"jobs"=>jobs);sorted=true),joinpath(out,"suite.toml"),"w")
    all(j["passed"] for j in jobs) || error("one or more jobs failed; inspect logs")
end
main()
