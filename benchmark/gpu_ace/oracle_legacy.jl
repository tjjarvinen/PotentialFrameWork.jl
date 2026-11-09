# Original fitted-model oracle in a separate dependency-compatible process.
using ACEpotentials, TOML, LinearAlgebra, SHA, BenchmarkTools
const J=ACEpotentials.JuLIP
const Old=ACEpotentials.ACE1
include("workloads.jl")
BLAS.set_num_threads(1)
source,destination=ARGS
ip=ACEpotentials.load_potential(source)
z=collect(ip.components[2].pibasis.zlist.list)
records=Dict[]
for g in (1,4)
    data=bcc_geometry(g,z)
    X=g==1 ? [.4 2.6 1.1;0. .2 2.8;0. .3 .1] : data.positions
    n=size(X,2); H=g==1 ? Matrix{Float64}(I,3,3)*16 : data.cell
    Z=g==1 ? [z[mod1(i,length(z))] for i=1:n] : data.species
    at=J.Atoms(X=[J.JVec(X[:,i]) for i=1:n],Z=Z,cell=H,pbc=(true,true,true))
    # SumIP methods are energy/forces/virial; each invokes original model code.
    e=J.energy(ip,at); f=reduce(hcat,J.forces(ip,at)); w=Matrix(J.virial(ip,at))
    oracle_call=()->(J.energy(ip,at),J.forces(ip,at),J.virial(ip,at))
    oracle_call()
    trial=run(@benchmarkable($oracle_call(),evals=1);samples=5,seconds=.2)
    push!(records,Dict("median_s"=>BenchmarkTools.median(trial).time*1e-9,"raw_ns"=>trial.times,"timed_boundary"=>"Three separate original JuLIP E/F/V calls; geometry and model resident","positions"=>vec(X),"cell"=>vec(H),"species"=>Int.(Z),"energy"=>e,"forces"=>vec(f),"virial"=>vec(w)))
end
open(destination,"w") do io
    TOML.print(io,Dict("source_sha256"=>bytes2hex(sha256(read(source))),"julia"=>string(VERSION),"threads"=>Threads.nthreads(),"blas_threads"=>BLAS.get_num_threads(),
        "ACEpotentials"=>string(pkgversion(ACEpotentials)),"records"=>records);sorted=true)
end
