# Dedicated pinned producer; no fitting, and never writes the source recipe.
# julia --project=data/ace_ref data/ace_ref/export_resolved_model.jl
using ACEpotentials, SHA, TOML, LinearAlgebra, SparseArrays
const M = ACEpotentials.Models
VERSION == v"1.13.2" || error("Use Julia 1.13.2")
pkgversion(ACEpotentials) == v"0.10.2" || error("Unsupported producer version")

include(joinpath(@__DIR__,"resolved_export.jl"))

path=joinpath(@__DIR__,"reference_model.json")
hash=bytes2hex(sha256(read(path)))
p,saved=ACEpotentials.load_model(path); model=p.model
isempty(p.ps.rbasis) && isempty(p.ps.pairbasis) || error("Unexpected trainable radial parameters")
flat=vcat(vec(p.ps.WB),vec(p.ps.Wpair))
flat==saved["model_parameters"] || error("Raw parameter ordering mismatch")
length(flat)==270 || error("Unexpected parameter count")
model.tensor.LL==(0,) || error("Only scalar invariants supported")
model.tensor.aabasis.hasconst && error("Constant products need a new schema")
model.Vref isa M.OneBody || error("Only one-body fixed references supported")
C=model.tensor.A2Bmaps[1]; rows,cols,vals=findnz(C)
versions=Dict(name=>saved["manifest"]["deps"][name][1]["version"] for name in
    ("ACEpotentials","EquivariantTensors","SpheriCart","Polynomials4ML","Lux","ACEfit"))
resolved=Dict("schema"=>"PotentialFrameWork.resolved_ACE", "schema_version"=>1,
    "source_sha256"=>hash,"producer_julia"=>string(VERSION),"producer_versions"=>versions,
    "species_atomic_numbers"=>collect(model._i2z),
    "reference_eV"=>[model.Vref.E0[z] for z in model._i2z],
    "angular"=>Dict("kind"=>"real_spherical", "normalization"=>"L2", "ordering"=>"l*l+l+m+1",
        "lmax"=>isqrt(length(model.ybasis))-1),
    "radial"=>radial_export(model.rbasis),"pair_radial"=>radial_export(model.pairbasis),
    "density_indices"=>[collect(t) for t in model.tensor.abasis.spec],
    "products"=>[collect(t) for spec in model.tensor.aabasis.specs for t in spec],
    "coupling"=>Dict("shape"=>collect(size(C)),"rows"=>rows,"columns"=>cols,"values"=>vals),
    "parameters"=>Dict("WB_shape"=>collect(size(p.ps.WB)),"Wpair_shape"=>collect(size(p.ps.Wpair)),
        "flat"=>flat,"ordering"=>"vec(WB), vec(Wpair); Julia column-major"))
open(joinpath(@__DIR__,"resolved_model_v1.toml"),"w") do io
    TOML.print(io,resolved;sorted=true)
end
bytes2hex(sha256(read(path)))==hash || error("Source model changed")
println("Exported 270 parameters, ",length(model.tensor.abasis)," densities, ",size(C)," coupling map")
