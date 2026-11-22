# Matched Float64, fresh-geometry public E/F/V on the supplied fit.
# Run separately in the package and reference environments.
using AtomsBase, AtomsCalculators, Unitful, LinearAlgebra, Statistics, TOML, SHA
include("../gpu_ace/workloads.jl")
const ROOT = normpath(joinpath(@__DIR__,"../.."))
const MODE,OUTPUT = ARGS
const DATA = joinpath(ROOT,"data/ace_ref")
const SOURCE = joinpath(DATA,"reference_model.json")
const SOURCE_SHA = bytes2hex(sha256(read(SOURCE)))
const SAMPLES = 9
BLAS.set_num_threads(1)

if MODE == "reference"
    using ACEpotentials
    calculator, _ = ACEpotentials.load_model(SOURCE)
    version = string(pkgversion(ACEpotentials))
elseif MODE == "pfw"
    using PotentialFrameWork, SpheriCart
    const A = PotentialFrameWork.ACE
    m, ps, st = A.load_resolved(joinpath(DATA,"resolved_model_v1.toml"))
    @assert m.source_sha256 == SOURCE_SHA
    calculator = A.ACECalculator(m,ps,st)
    version = string(pkgversion(PotentialFrameWork))
else
    error("expected reference or pfw")
end

d = bcc_geometry(2,[13,22])
frames = [AtomsBase.atomic_system(
    [ChemicalSpecies(d.species[i]) => (d.positions[:,i] .+ a .* sin.(i .+ (1:3))) .* u"Å" for i in 1:16],
    ntuple(j -> d.cell[:,j]*u"Å",3), (true,true,true)) for a in (0.,1e-5,-1e-5)]
run(sys) = MODE == "reference" ? AtomsCalculators.energy_forces_virial(sys,calculator) :
    AtomsCalculators.energy_forces_virial(sys,calculator;executor=:threaded,workers=Threads.nthreads(:default))
numeric(r) = Dict("energy"=>ustrip(u"eV",r.energy),
    "forces"=>[collect(ustrip.(u"eV/Å",x)) for x in r.forces],
    "virial"=>[collect(ustrip.(u"eV",r.virial[:,j])) for j in 1:3])
try
    expected = numeric.(run.(frames))
    foreach(run, frames)
    GC.gc()
    trials = [@timed run(frames[mod1(i,3)]) for i in 1:SAMPLES]
    times = getproperty.(trials,:time)
    report = Dict("mode"=>MODE,"julia"=>string(VERSION),"version"=>version,
        "public_boundary"=>"fresh-geometry Float64 E/F/V with owning outputs",
        "source_sha256"=>SOURCE_SHA,"harness_sha256"=>bytes2hex(sha256(read(@__FILE__))),
        "cpu"=>Sys.cpu_info()[1].model,"threads"=>Threads.nthreads(:default),
        "affinity"=>strip(read(`taskset -pc $(getpid())`,String)),"blas_threads"=>BLAS.get_num_threads(),
        "atoms"=>16,"seconds"=>times,"median_s"=>median(times),"min_s"=>minimum(times),"max_s"=>maximum(times),
        "bytes"=>getproperty.(trials,:bytes),"numeric_outputs"=>expected)
    open(io -> TOML.print(io,report;sorted=true),OUTPUT,"w")
finally
    MODE == "pfw" && close(calculator)
end
