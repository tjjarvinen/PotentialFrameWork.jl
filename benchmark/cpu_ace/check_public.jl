# Separate correctness gate for the benchmark's new model/size/precision grid.
include("common.jl")
function main()
    data,oracles,output=ARGS; records=Dict[]
    for name in ("small-tial","medium-ti","large-tial","ti","tialv")
        m,ps,_=A.load_resolved(joinpath(data,name*"-resolved.toml"))
        oracle=TOML.parsefile(joinpath(oracles,name*"-oracle.toml"))
        for g in (name in ("ti","tialv") ? (4,) : (4,8,18)), T in (Float64,Float32)
            # Large Float32 model at 10k is checked numerically, not timed here.
            b=A.adapt_model(m,ps;precision=T); sys=structure(bcc_geometry(g,m.atomic_numbers);T)
            raw=only(filter(r->r["g"]==g,oracle["records"]))
            ref=(;energy=raw["energy"],forces=reshape(raw["forces"],3,:),virial=reshape(raw["virial"],3,3))
            for algorithm in (:reference,:optimized)
                calc=A.ACECalculator(b.model,b.parameters;precision=T,product_algorithm=algorithm)
                public=AC.energy_forces_virial(sys,calc)
                y=(;energy=U.ustrip(U.u"eV",public.energy),
                    forces=reduce(hcat,[U.ustrip.(U.u"eV/Å",v) for v in public.forces]),
                    virial=Matrix(U.ustrip.(U.u"eV",public.virial)))
                agree(y,ref,T)
                numeric=A.prepared_result(calc.executor); agree(y,numeric,T)
                push!(records,Dict("model"=>name,"atoms"=>length(sys),"precision"=>string(T),"algorithm"=>string(algorithm),
                    "oracle_version"=>oracle["ACEpotentials"],"external_errors"=>accuracy(y,ref),"prepared_errors"=>accuracy(y,numeric)))
                open(io->TOML.print(io,Dict("schema"=>"PFW.CPUACE.PublicValidation.v1","julia"=>string(VERSION),"cases"=>records);sorted=true),output,"w")
            end
        end
    end
    println("Public E/F/V checks passed for ",length(records)," model/size/precision/algorithm cases")
end
main()
