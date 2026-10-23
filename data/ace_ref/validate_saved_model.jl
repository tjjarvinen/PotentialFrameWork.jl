# Run with: julia --project=data/ace_ref data/ace_ref/validate_saved_model.jl
# Deliberately does not include generate_reference_data.jl or fit a model.
using ACEpotentials, ExtXYZ, AtomsCalculators, Unitful
using SHA, TOML, LinearAlgebra

const fixture_dir = @__DIR__
const model_path = joinpath(fixture_dir, "reference_model.json")
const xyz_path = joinpath(fixture_dir, "reference_results.xyz")
const model_hash = bytes2hex(sha256(read(model_path)))
VERSION == v"1.13.2" || error("Use the pinned Julia 1.13.2 reference environment")
pkgversion(ACEpotentials) == v"0.10.2" || error("Unexpected ACEpotentials version")
model, saved = ACEpotentials.load_model(model_path)
frames = ExtXYZ.load(xyz_path)
length(frames) == 33 || error("Expected all 33 reference frames")
errors = zeros(3)
records = map(enumerate(frames)) do (i, sys)
    xyz = ExtXYZ.Atoms(sys)
    haskey(xyz.atom_data, :force) || error("Frame $i lacks per-atom prediction forces")
    result = AtomsCalculators.energy_forces_virial(sys, model)
    energy = ustrip(u"eV", result.energy)
    forces = [collect(ustrip.(u"eV/Å", f)) for f in result.forces]
    virial = ustrip.(u"eV", result.virial)
    errors[1] = max(errors[1], abs(energy - xyz.system_data.energy))
    errors[2] = max(errors[2], maximum(abs, hcat(forces...) - hcat(xyz.atom_data.force...)))
    errors[3] = max(errors[3], maximum(abs, virial - xyz.system_data.virial))
    Dict("frame" => i, "energy_eV" => energy, "forces_eV_per_A" => forces,
         "virial_eV_columns" => [collect(virial[:,a]) for a=1:3])
end
# These bounds describe the supplied eight-decimal geometry/output fixture.
# Full precision comparisons must use the recomputed records, not these bounds.
all(errors .< [1e-5, 1e-5, 1e-4]) || error("Rounded fixture mismatch: $errors")
bytes2hex(sha256(read(model_path))) == model_hash || error("Saved model was modified")
output = Dict("schema_version" => 1, "model_sha256" => model_hash,
    "geometry_sha256" => bytes2hex(sha256(read(xyz_path))),
    "julia_version" => string(VERSION), "ACEpotentials_version" => string(pkgversion(ACEpotentials)),
    "geometry" => "reloaded reference_results.xyz; original atom order",
    "maximum_rounded_errors_E_F_W" => errors, "frames" => records)
open(joinpath(fixture_dir, "reloaded_predictions.toml"), "w") do io
    TOML.print(io, output)
end
println("Validated ", length(frames), " frames; max rounded E/F/W differences = ", errors)
println("model SHA256 = ", model_hash)
