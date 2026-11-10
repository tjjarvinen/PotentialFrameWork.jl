# Original analytic ACE1, never labelled modern. Run in the legacy environment.
using ACEpotentials, TOML, LinearAlgebra, SHA, BenchmarkTools
const J=ACEpotentials.JuLIP
pkgversion(ACEpotentials)==v"0.6.6" || error("requires isolated ACEpotentials 0.6.6")
include(joinpath(@__DIR__,"../gpu_ace/workloads.jl"))
BLAS.set_num_threads(1)
source,destination=ARGS
ip=ACEpotentials.load_potential(source)
z=collect(ip.components[2].pibasis.zlist.list)
sizes=parse.(Int,split(get(ENV,"PFW_CPU_SIZES","4"),','))
samples=parse(Int,get(ENV,"PFW_CPU_SAMPLES","10")); seconds=parse(Float64,get(ENV,"PFW_CPU_SECONDS","3"))
records=Dict[]
for g in sizes
    d=bcc_geometry(g,z); n=size(d.positions,2)
    at=J.Atoms(X=[J.JVec(d.positions[:,i]) for i=1:n],Z=d.species,cell=d.cell,pbc=(true,true,true))
    f=()->(J.energy(ip,at),reduce(hcat,J.forces(ip,at)),Matrix(J.virial(ip,at)))
    first=@elapsed r=f()
    @assert isfinite(r[1]) && all(isfinite,r[2]) && all(isfinite,r[3])
    trial=run(@benchmarkable($f(),evals=1);samples,seconds)
    push!(records,Dict("g"=>g,"atoms"=>n,"positions"=>vec(d.positions),"cell"=>vec(d.cell),"species"=>Int.(d.species),
        "energy"=>r[1],"forces"=>vec(r[2]),"virial"=>vec(r[3]),"first_call_s"=>first,
        "raw_ns"=>trial.times,"host_bytes"=>trial.memory,"samples"=>length(trial.times),"undersampled"=>length(trial.times)<10,
        "boundary"=>"legacy ACEpotentials 0.6.6 / ACE1 analytic; three separate JuLIP E/F/V calls; native fixed-geometry neighbor cache; owning outputs"))
end
open(destination,"w") do io
    TOML.print(io,Dict("schema"=>"PFW.CPUACE.LegacyOracle.v1","source_sha256"=>bytes2hex(sha256(read(source))),
        "julia"=>string(VERSION),"ACEpotentials"=>string(pkgversion(ACEpotentials)),"threads"=>Threads.nthreads(),"blas_threads"=>BLAS.get_num_threads(),"records"=>records);sorted=true)
end
