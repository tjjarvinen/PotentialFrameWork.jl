# Isolated modern ACEpotentials producer. No fitting and no existing input writes.
using ACEpotentials, SHA, TOML, LinearAlgebra, SparseArrays, BenchmarkTools
using AtomsBase, AtomsCalculators, Unitful, StaticArrays
const M=ACEpotentials.Models
pkgversion(ACEpotentials)==v"0.10.2" || error("requires modern ACEpotentials 0.10.2")
# Reuse the independently validated exact spline exporter, without its script body.
exporter=joinpath(@__DIR__,"../../data/ace_ref/export_resolved_model.jl")
src=read(exporter,String)
include_string(@__MODULE__,split(split(src,"function plain_fields";limit=2)[2],"\npath=joinpath";limit=2)[1] |> s->"function plain_fields"*s,exporter)
include(joinpath(@__DIR__,"../gpu_ace/workloads.jl"))
BLAS.set_num_threads(1)
out=abspath(ARGS[1]); mkpath(out)
sizes=parse.(Int,split(get(ENV,"PFW_CPU_SIZES","4,8,18"),','))
samples=parse(Int,get(ENV,"PFW_CPU_SAMPLES","10")); seconds=parse(Float64,get(ENV,"PFW_CPU_SECONDS","3"))
for (name,elements,order,degree) in (("medium-ti",[:Ti],3,10),("large-tial",[:Ti,:Al],4,12),("small-tial",[:Ti,:Al],3,6))
    println("Generating ",name); flush(stdout)
    p= name=="small-tial" ? first(ACEpotentials.load_model(joinpath(@__DIR__,"../../data/ace_ref/reference_model.json"))) :
        ACEpotentials.ace1_model(elements=elements,order=order,totaldegree=degree,rcut=5.5,Eref=Dict(z=>0. for z in elements))
    model=p.model
    if name!="small-tial"
        # Every central species has deterministic nonzero, decorrelated weights.
        for z in axes(p.ps.WB,2), i in axes(p.ps.WB,1)
            p.ps.WB[i,z]=.02*sin(i*1.7+z*.91)/(1+i/40)
        end
        for z in axes(p.ps.Wpair,2),i in axes(p.ps.Wpair,1)
            p.ps.Wpair[i,z]=.01*cos(i*1.3+z*.7)/(1+i)
        end
    end
    recipe=joinpath(out,name*".json")
    ACEpotentials.save_model(p,recipe;verbose=false,meta=Dict("synthetic"=>name!="small-tial","physical_accuracy_claim"=>false))
    C=model.tensor.A2Bmaps[1]; rows,cols,vals=findnz(C)
    resolved=Dict("schema"=>"PotentialFrameWork.resolved_ACE","schema_version"=>1,
        "source_sha256"=>bytes2hex(sha256(read(recipe))),"producer_julia"=>string(VERSION),
        "producer_versions"=>Dict("ACEpotentials"=>string(pkgversion(ACEpotentials))),
        "species_atomic_numbers"=>collect(model._i2z),"reference_eV"=>[model.Vref.E0[z] for z in model._i2z],
        "angular"=>Dict("kind"=>"real_spherical","normalization"=>"L2","ordering"=>"l*l+l+m+1","lmax"=>isqrt(length(model.ybasis))-1),
        "radial"=>radial_export(model.rbasis),"pair_radial"=>radial_export(model.pairbasis),
        "density_indices"=>[collect(t) for t in model.tensor.abasis.spec],
        "products"=>[collect(t) for spec in model.tensor.aabasis.specs for t in spec],
        "coupling"=>Dict("shape"=>collect(size(C)),"rows"=>rows,"columns"=>cols,"values"=>vals),
        "parameters"=>Dict("WB_shape"=>collect(size(p.ps.WB)),"Wpair_shape"=>collect(size(p.ps.Wpair)),
            "flat"=>vcat(vec(p.ps.WB),vec(p.ps.Wpair)),"ordering"=>"vec(WB), vec(Wpair); Julia column-major"))
    open(joinpath(out,name*"-resolved.toml"),"w") do io; TOML.print(io,resolved;sorted=true); end
    records=Dict[]
    for g in sizes
        d=bcc_geometry(g,collect(model._i2z))
        sys=AtomsBase.atomic_system([ChemicalSpecies(d.species[i])=>SVector{3}(d.positions[:,i])*u"Å" for i=1:size(d.positions,2)],ntuple(j->SVector{3}(d.cell[:,j])*u"Å",3),(true,true,true))
        f=()->AtomsCalculators.energy_forces_virial(sys,p)
        elapsed=@elapsed r=f()
        e=ustrip(u"eV",r.energy); F=reduce(hcat,[ustrip.(u"eV/Å",v) for v in r.forces]); W=Matrix(ustrip.(u"eV",r.virial))
        all(isfinite,F) && all(isfinite,W) && isfinite(e) || error("nonfinite fixture")
        maximum(abs,F)<1e4 && abs(e)/length(sys)<1e4 || error("ill-scaled fixture")
        trial=run(@benchmarkable($f(),evals=1);samples,seconds)
        push!(records,Dict("g"=>g,"atoms"=>length(sys),"positions"=>vec(d.positions),"cell"=>vec(d.cell),"species"=>d.species,
            "energy"=>e,"forces"=>vec(F),"virial"=>vec(W),"raw_ns"=>trial.times,"host_bytes"=>trial.memory,
            "first_call_s"=>elapsed,"samples"=>length(trial.times),"undersampled"=>length(trial.times)<10,
            "boundary"=>"modern ACEpotentials public combined E/F/V; fresh neighbors; owning outputs; Float64"))
        println(name," ",length(sys)," atoms ",median(trial).time/1e6," ms"); flush(stdout)
    end
    report=Dict("schema"=>"PFW.CPUACE.ModernOracle.v1","model"=>name,"synthetic"=>name!="small-tial",
        "backend_type"=>string(typeof(model)),"tensor_type"=>string(typeof(model.tensor)),
        "julia"=>string(VERSION),"ACEpotentials"=>string(pkgversion(ACEpotentials)),"threads"=>Threads.nthreads(),"blas_threads"=>BLAS.get_num_threads(),
        "source_sha256"=>resolved["source_sha256"],"resolved_sha256"=>bytes2hex(sha256(read(joinpath(out,name*"-resolved.toml")))),
        "precision_unsupported"=>["Float32: producer/calculator uses Float64 tables; no matching precision claimed"],"records"=>records)
    open(joinpath(out,name*"-oracle.toml"),"w") do io; TOML.print(io,report;sorted=true); end
end
