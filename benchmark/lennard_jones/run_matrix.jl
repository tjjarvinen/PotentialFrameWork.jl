# Run from the root with --project=benchmark/lennard_jones after setup.
# Workers run sequentially, so benchmarks at different thread counts do not compete.
counts=parse.(Int,split(get(ENV,"PFW_LJ_THREADS","1,2,4,8"),','))
all(>(0),counts) || error("thread counts must be positive")
output=isempty(ARGS) ? joinpath(@__DIR__,"results") : abspath(ARGS[1])
mkpath(output)
for n in unique(counts)
    println("Starting independent Julia process with $n threads"); flush(stdout)
    cmd=`$(Base.julia_cmd()) --project=$(@__DIR__) --threads=$n $(joinpath(@__DIR__,"worker.jl")) $output`
    run(cmd)
end
println("CSV results and exact dependency versions: $output")
