# Existing independently reconstructed ACEpotentials 0.10.2 model; never refit.
using LinearAlgebra, ACEpotentials, AtomsBase, AtomsCalculators, Unitful, TOML, SHA, BenchmarkTools, StaticArrays
include("workloads.jl")
BLAS.set_num_threads(1)
source,destination=ARGS
pkgversion(ACEpotentials)==v"0.10.2" || error("use the supplied pinned TiAl reference environment")
oracle,_=ACEpotentials.load_model(source)
records=Dict[]
for g in (1,4)
    data=bcc_geometry(g,[22,13])
    X=g==1 ? [.4 2.6 1.1;0. .2 2.8;0. .3 .1] : data.positions
    H=g==1 ? Matrix{Float64}(I,3,3)*16 : data.cell
    n=size(X,2); Z=[[22,13][mod1(i,2)] for i=1:n]
    sys=AtomsBase.atomic_system([ChemicalSpecies(Z[i])=>SVector{3}(X[:,i])*u"Å" for i=1:n],ntuple(j->SVector{3}(H[:,j])*u"Å",3),(true,true,true))
    f=()->AtomsCalculators.energy_forces_virial(sys,oracle)
    r=f(); trial=run(@benchmarkable($f(),evals=1);samples=5,seconds=.2)
    push!(records,Dict("positions"=>vec(X),"cell"=>vec(H),"species"=>Z,"energy"=>ustrip(u"eV",r.energy),
        "forces"=>vec(reduce(hcat,[ustrip.(u"eV/Å",v) for v in r.forces])),"virial"=>collect(vec(ustrip.(u"eV",r.virial))),
        "median_s"=>BenchmarkTools.median(trial).time*1e-9,"raw_ns"=>trial.times,"timed_boundary"=>"Original ACEpotentials AtomsCalculators energy_forces_virial public call"))
end
open(destination,"w") do io; TOML.print(io,Dict("source_sha256"=>bytes2hex(sha256(read(source))),"julia"=>string(VERSION),"threads"=>Threads.nthreads(),"blas_threads"=>BLAS.get_num_threads(),"ACEpotentials"=>string(pkgversion(ACEpotentials)),"records"=>records);sorted=true); end
