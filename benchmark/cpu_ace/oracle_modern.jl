# Re-time cached recipes without rewriting model/source/export artifacts.
using ACEpotentials, AtomsBase, AtomsCalculators, Unitful, StaticArrays
using BenchmarkTools, LinearAlgebra, TOML, SHA
include(joinpath(@__DIR__,"../gpu_ace/workloads.jl"))
function main()
    pkgversion(ACEpotentials)==v"0.10.2" || error("requires modern ACEpotentials 0.10.2")
    data,out=abspath.(ARGS); mkpath(out); BLAS.set_num_threads(1)
    sizes=parse.(Int,split(get(ENV,"PFW_CPU_SIZES","4,8,18"),','))
    samples=parse(Int,get(ENV,"PFW_CPU_SAMPLES","10")); seconds=parse(Float64,get(ENV,"PFW_CPU_SECONDS","60"))
    for name in ("small-tial","medium-ti","large-tial")
        recipe=joinpath(data,name*".json"); p,_=ACEpotentials.load_model(recipe); records=Dict[]
        for g in sizes
            d=bcc_geometry(g,collect(p.model._i2z)); n=size(d.positions,2)
            sys=AtomsBase.atomic_system([ChemicalSpecies(d.species[i])=>SVector{3}(d.positions[:,i])*u"Å" for i=1:n],ntuple(j->SVector{3}(d.cell[:,j])*u"Å",3),(true,true,true))
            f=()->AtomsCalculators.energy_forces_virial(sys,p)
            first=@elapsed r=f()
            F=reduce(hcat,[ustrip.(u"eV/Å",v) for v in r.forces]); W=Matrix(ustrip.(u"eV",r.virial)); e=ustrip(u"eV",r.energy)
            @assert all(isfinite,F) && all(isfinite,W) && isfinite(e)
            trial=run(@benchmarkable($f(),evals=1);samples,seconds)
            push!(records,Dict("g"=>g,"atoms"=>n,"positions"=>vec(d.positions),"cell"=>vec(d.cell),"species"=>d.species,
                "energy"=>e,"forces"=>vec(F),"virial"=>vec(W),"first_call_s"=>first,
                "raw_ns"=>trial.times,"host_bytes"=>trial.memory,"samples"=>length(trial.times),"undersampled"=>length(trial.times)<10,
                "boundary"=>"modern ACEpotentials 0.10.2 public combined E/F/V; fresh neighbors; owning outputs; Float64"))
            println(name," ",n," ",median(trial).time/1e6," ms; ",length(trial.times)," samples"); flush(stdout)
        end
        report=Dict("schema"=>"PFW.CPUACE.ModernOracle.v1","model"=>name,"source_sha256"=>bytes2hex(sha256(read(recipe))),
            "backend_type"=>string(typeof(p.model)),"tensor_type"=>string(typeof(p.model.tensor)),
            "julia"=>string(VERSION),"ACEpotentials"=>string(pkgversion(ACEpotentials)),"threads"=>Threads.nthreads(),"blas_threads"=>BLAS.get_num_threads(),
            "precision_unsupported"=>["Float32: Float64 radial/model tables; no external Float32 comparison claimed"],"records"=>records)
        open(io->TOML.print(io,report;sorted=true),joinpath(out,name*"-oracle.toml"),"w")
    end
end
main()
