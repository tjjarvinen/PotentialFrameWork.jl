# Local setup, in a fresh Julia environment:
# julia --project=. -e 'using Pkg; Pkg.add(Pkg.PackageSpec(name="ACEpotentials", version=v"0.10.2")); Pkg.pin("ACEpotentials"); Pkg.add(["ACEfit", "AtomsBase", "ExtXYZ", "Lux", "Unitful"])'
#
# Inspect the full fit without assembling its design matrix:
# julia --project=. -t auto train_may30_ace.jl ti64_combined_data_May30.extxyz
#
# Smoke fit: degree 8, 200 training structures, QR:
# julia --project=. -t auto train_may30_ace.jl ti64_combined_data_May30.extxyz smoke.jls 8 200 qr
#
# Full fit: degree 15, all non-held-out structures, BLR:
# julia --project=. -t auto train_may30_ace.jl ti64_combined_data_May30.extxyz may30.jls
#
# Reload in the same Julia/package environment:
# include("train_may30_ace.jl")
# calc, checkpoint = load_may30_fit("may30.jls")
# etace = ACEpotentials.ETModels.convert2et_full(calc.model, calc.ps, calc.st)

using ACEpotentials, ACEfit, AtomsBase, ExtXYZ, Lux, Unitful
using LinearAlgebra, Random, Serialization, SHA

const M = ACEpotentials.Models
const ELEMENTS = (:Ti, :Al, :V)
const ATOMIC_NUMBERS = (22, 13, 23)
const SEED = 1234
const LABELS = (energy_key="QM_energy", force_key="QM_forces", virial_key="QM_virial")

function build_model(settings)
    cutoffs = map(M._default_rin0cuts(ELEMENTS)) do entry
        merge(entry, (rcut=settings.rcut,))
    end
    e0s = Dict(element => value * u"eV" for (element, value) in zip(ELEMENTS, settings.e0))
    model = M.ace_model(;
        elements=ELEMENTS, order=settings.order, max_level=settings.degree,
        level=M.TotalDegree(Float64(length(ELEMENTS)), 1 / 1.5),
        pair_maxn=settings.degree, rin0cuts=cutoffs, E0s=e0s,
        init_Wradial=:onehot, pair_learnable=true,
    )
    ps, st = Lux.setup(MersenneTwister(SEED), model)
    return ACEpotentials.ACEPotential(model, ps, st)
end

function check_version()
    pkgversion(ACEpotentials) == v"0.10.3" || error("Use ACEpotentials 0.10.3")
end

function load_may30_fit(path)
    check_version()
    saved = deserialize(path)
    saved.schema == 1 || error("Unsupported May30 checkpoint schema")
    model = build_model(saved.settings).model
    return ACEpotentials.ACEPotential(model, saved.parameters, saved.state), saved
end

function labels_for(data)
    labels = map(enumerate(data)) do (i, sys)
        record = ACEpotentials.AtomsData(sys; LABELS...)
        all(!isnothing, (record.energy_key, record.force_key, record.virial_key)) ||
            error("Frame $i is missing a QM energy, forces or virial field")
        record
    end
    return labels
end

function reference_energies(records, train_indices)
    composition = zeros(length(train_indices), length(ELEMENTS))
    energies = zeros(length(train_indices))
    for (row, index) in enumerate(train_indices)
        record = records[index]
        sys = record.system
        for atom in 1:length(sys)
            z = Int(AtomsBase.species(sys, atom).atomic_number)
            col = findfirst(==(z), ATOMIC_NUMBERS)
            isnothing(col) && error("Unexpected atomic number $z in frame $index")
            composition[row, col] += 1
        end
        composition[row, :] ./= length(sys)
        energies[row] = ustrip(sys[record.energy_key]) / length(sys)
    end
    rank(composition) == length(ELEMENTS) ||
        error("Training compositions cannot determine three reference energies. Use more frames.")
    return Tuple(composition \ energies)
end

function main(args)
    check_version()
    1 <= length(args) <= 5 ||
        error("Usage: train_may30_ace.jl DATA.extxyz [OUTPUT.jls [DEGREE [MAX_TRAIN [blr|qr]]]]")
    dataset = args[1]
    output = get(args, 2, nothing)
    degree = parse(Int, get(args, 3, "15"))
    max_train = parse(Int, get(args, 4, string(typemax(Int))))
    solver_name = get(args, 5, "blr")
    degree > 0 && max_train > 0 || error("DEGREE and MAX_TRAIN must be positive")
    solver_name in ("blr", "qr") || error("Solver must be blr or qr")
    if !isnothing(output)
        ispath(output) && error("Output already exists: $output")
        isdir(dirname(abspath(output))) || error("Output directory does not exist")
    end

    data = ExtXYZ.load(dataset)
    length(data) >= 10 || error("Expected at least ten configurations")
    records = labels_for(data)
    shuffled = randperm(MersenneTwister(SEED), length(data))
    nheld = max(1, length(data) ÷ 10)
    heldout_indices = sort(shuffled[1:nheld])
    ntrain = min(max_train, length(data) - nheld)
    train_indices = sort(shuffled[nheld+1:nheld+ntrain])
    e0 = reference_energies(records, train_indices)

    settings = (order=3, degree=degree, rcut=6.0, e0=e0)
    calc = build_model(settings)
    rows = sum(3 * length(data[i]) + 7 for i in train_indices)
    columns = M.length_basis(calc.model)
    gib = Float64(rows) * columns * sizeof(Float64) / 2.0^30
    println("Frames: $ntrain train, $nheld held out (seed $SEED)")
    println("Order 3, degree $degree, cutoff 6 Å. Reference energies (Ti, Al, V): $e0 eV")
    println("Fit matrix: $rows × $columns, approximately $(round(gib; digits=2)) GiB per dense Float64 copy.")
    println("Assembly, weighting and the solver can require multiple copies plus workspace.")
    isnothing(output) && return

    weights = Dict("default" => Dict("E"=>30.0, "F"=>1.0, "V"=>1.0))
    solver = solver_name == "blr" ? ACEfit.BLR() : ACEfit.QR(lambda=1e-3)
    ACEpotentials.acefit!(data[train_indices], calc; LABELS..., weights, solver, smoothness=4)
    all(isfinite, calc.ps.WB) && all(isfinite, calc.ps.Wpair) || error("Fit produced nonfinite coefficients")
    checkpoint = (schema=1, settings=settings, parameters=calc.ps, state=calc.st,
                  train_indices=train_indices, heldout_indices=heldout_indices,
                  solver=solver_name, weights=weights, smoothness=4,
                  data_sha256=open(io -> bytes2hex(sha256(io)), dataset),
                  julia_version=VERSION, acepotentials_version=pkgversion(ACEpotentials))
    destination = abspath(output)
    mktemp(dirname(destination)) do temporary, io
        serialize(io, checkpoint)
        close(io)
        mv(temporary, destination)
    end
    println("Saved fitted polynomial checkpoint: $destination")
    println("Held-out frame errors (random split for conversion/benchmark development):")
    ACEpotentials.compute_errors(data[heldout_indices], calc; LABELS..., weights)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main(ARGS)
end
