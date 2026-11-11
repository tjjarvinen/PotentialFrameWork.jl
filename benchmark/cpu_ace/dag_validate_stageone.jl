# Export exact full-model outputs from immutable stage one, then compare current
# polynomial/prefix/DAG on identical perturbed fixtures. No timing claims here.
include("common.jl")
function main()
    mode,path=ARGS; BLAS.set_num_threads(1)
    if mode=="export"
        records=Dict[]
        for name in ("small-tial","medium-ti","large-tial","ti","tialv"), g in (4,8,18), T in (Float64,Float32)
            g==18 && !(name in ("small-tial","large-tial")) && continue
            m,ps,_=A.load_resolved(joinpath(@__DIR__,"results/2026-11-10/models",name*"-resolved.toml");precision=T)
            sys=structure(bcc_geometry(g,m.atomic_numbers);T)
            input=A.ace_inputs(sys,m;environment=A.FullEnvironment(width=4),precision=T)
            w=A.CPUExecutor(m,input.environment;precision=T,product_algorithm=:optimized)
            A.execute_prepared!(w,input.positions,input.cell,ps); y=A.prepared_result(w)
            push!(records,Dict("model"=>name,"g"=>g,"precision"=>string(T),"geometry_sha256"=>bytes2hex(sha256(reinterpret(UInt8,vec(input.positions)))),
                "energy"=>y.energy,"forces"=>vec(y.forces),"virial"=>vec(y.virial)))
        end
        open(io->TOML.print(io,Dict("source_files_sha256"=>source_provenance(dirname(dirname(pathof(PotentialFrameWork)))),"records"=>records);sorted=true),path,"w")
    elseif mode=="compare"
        reference=TOML.parsefile(path); records=Dict[]
        for raw in reference["records"]
            name=raw["model"]; T=raw["precision"]=="Float32" ? Float32 : Float64
            m,ps,_=A.load_resolved(joinpath(@__DIR__,"results/2026-11-10/models",name*"-resolved.toml");precision=T)
            sys=structure(bcc_geometry(raw["g"],m.atomic_numbers);T)
            input=A.ace_inputs(sys,m;environment=A.FullEnvironment(width=4),precision=T)
            @assert raw["geometry_sha256"]==bytes2hex(sha256(reinterpret(UInt8,vec(input.positions))))
            w=A.CPUExecutor(m,input.environment;precision=T)
            ref=(;energy=raw["energy"],forces=reshape(raw["forces"],3,:),virial=reshape(raw["virial"],3,3))
            for algorithm in (:reference,:optimized,:dag)
                A.execute_prepared!(w,input.positions,input.cell,ps;product_algorithm=algorithm)
                y=A.prepared_result(w); agree(y,ref,T)
                push!(records,Dict("model"=>name,"atoms"=>length(sys),"precision"=>string(T),"algorithm"=>string(algorithm),"stageone_errors"=>accuracy(y,ref)))
            end
        end
        out=joinpath(dirname(path),"stageone-comparison.toml")
        open(io->TOML.print(io,Dict("source_files_sha256"=>source_provenance(dirname(dirname(pathof(PotentialFrameWork)))),"stageone_output_sha256"=>bytes2hex(sha256(read(path))),"cases"=>records);sorted=true),out,"w")
        println("Immutable stage-one full E/F/V comparisons passed: ",length(records))
    else
        error("mode: export or compare")
    end
end
main()
