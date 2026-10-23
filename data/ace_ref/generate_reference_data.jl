using ACEpotentials
using ExtXYZ
using AtomsCalculators
using Unitful

data, _, meta = ACEpotentials.example_dataset("TiAl_tutorial")
train_data = data[1:5:end];
test_data = data[2:10:end];

hyperparams = (elements = [:Ti, :Al],
					order = 3,
					totaldegree = 6,
					rcut = 5.5,
					Eref = [:Ti => -1586.0195, :Al => -105.5954])
model = ace1_model(; hyperparams...)

weights = Dict(
        "FLD_TiAl" => Dict("E" => 60.0, "F" => 1.0 , "V" => 1.0 ),
        "TiAl_T5000" => Dict("E" => 5.0, "F" => 1.0 , "V" => 1.0 ));

solver = ACEfit.BLR()

P = algebraic_smoothness_prior(model; p = 4)    #  (p = 4 is in fact the default)

ace_calc = acefit!(train_data, model; solver=solver, prior = P, weights=weights);

ACEpotentials.save_model(ace_calc, "reference_model.json")

refs = map(test_data) do sys
    result = AtomsCalculators.energy_forces_virial(sys, ace_calc)

    e = ustrip(u"eV", result.energy)
    f = [ustrip.(u"eV/Å", force) for force in result.forces]
    v = ustrip.(u"eV", result.virial)

    xyz = ExtXYZ.Atoms(sys)

    # Replace the dataset's per-atom forces with model predictions.
    atom_data = merge(xyz.atom_data, (; force=f))
    system_data = merge(xyz.system_data, (; energy=e, virial=v))

    ExtXYZ.Atoms(atom_data, system_data)
end


ExtXYZ.save("reference_results.xyz", refs)