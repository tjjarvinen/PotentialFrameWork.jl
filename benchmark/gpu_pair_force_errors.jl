# Replay deterministic full-force cases to add per-atom vector-norm errors to an
# existing timing dataset; preserve all timing fields. The main benchmark also
# records these fields on future runs.
prefix=first(split(read(joinpath(@__DIR__,"gpu_pairs.jl"),String),"\nfor m in (14,29)"))
include_string(@__MODULE__,prefix,joinpath(@__DIR__,"gpu_pairs.jl"))
path=joinpath(@__DIR__,"gpu_pair_measurements.toml")
data=TOML.parsefile(path)
cache=Dict{Tuple,Any}()
for row in data["measurements"]
    if row["request"]=="EnergyOnly"
        row["max_atom_force_error_eV_A"]=0.
        row["max_atom_force_error_original_eV_A"]=0.
        continue
    end
    T=row["precision"]=="Float32" ? Float32 : Float64
    key=(row["replication"],Symbol(row["kind"]),T)
    sys,p,ref,original=get!(cache,key) do
        m,kind,_=key; sys=structure(m,kind)
        p=kind===:argon ? PP.lennard_jones() : SmoothPair([:Ar,:Ne],[.0103 .008;.008 .006],8.5)
        (sys,p,evaluate(sys,p;precision=T,parallel=true),evaluate(sys,p;precision=Float64,parallel=true))
    end
    y=PP.evaluate_gpu(sys,p,backend;precision=T,config=NL.GPUConfig(;strategy=Symbol(row["strategy"])))
    delta=ustrip.(y.forces-ref.forces)
    row["max_atom_force_error_eV_A"]=maximum(sqrt.(sum(abs2,delta;dims=1)))
    delta=ustrip.(y.forces)-ustrip.(original.forces)
    row["max_atom_force_error_original_eV_A"]=maximum(sqrt.(sum(abs2,delta;dims=1)))
end
data["force_norm_validation"]="separate deterministic replay; timing fields unchanged"
open(io->TOML.print(io,data),path,"w")
