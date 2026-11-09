function source_provenance(root)
    paths=String[joinpath(root,"Project.toml")]
    for dir in ("src","ext")
        for (base,_,files) in walkdir(joinpath(root,dir)), file in files
            endswith(file,".jl") && push!(paths,joinpath(base,file))
        end
    end
    for file in ("run.jl","workloads.jl","audit.jl","suite.jl")
        push!(paths,joinpath(@__DIR__,file))
    end
    Dict(relpath(path,root)=>bytes2hex(SHA.sha256(read(path))) for path in sort(paths))
end
