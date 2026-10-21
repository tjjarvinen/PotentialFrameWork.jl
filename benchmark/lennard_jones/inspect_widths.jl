# Inspect both the callback and the fused cell traversal on the executing CPU.
using PotentialFrameWork, AtomsBase, InteractiveUtils
include("common.jl")
using .ArgonComparison
const NL=PotentialFrameWork.NeighbourList
const PP=PotentialFrameWork.PairPotentials

function inspect_widths()
    directory=isempty(ARGS) ? joinpath(@__DIR__,"codegen") : abspath(ARGS[1])
    mkpath(directory)
    p,_=argon_models()
    ar=ChemicalSpecies(:Ar)
    summary=String["# LJ vector-width code inspection", "",
        "Julia $(VERSION), CPU $(Sys.cpu_info()[1].model), target $(Sys.CPU_NAME), Float64.", "",
        "Generated on the executing CPU; inspect the native files before drawing", 
        "conclusions about another CPU target. W=8 carries eight Float64 lanes (512 bits).", "",
        "| Width | Target | LLVM vector multiply lines | ZMM packed-double arithmetic lines | YMM packed-double arithmetic lines |",
        "|---:|---|---:|---:|---:|"]
    for width in (1,4,8)
        w=Workspace(;width)
        pack(f)=NL.pack(Val(width),f)
        block=PairBlock(ar,ntuple(_ -> ar,width),pack(k -> 3.7+0.01*k),
                       pack(_ -> 0.2),pack(_ -> 0.1),pack(_ -> true),pack(_ -> true))
        targets=(("callback",p,Tuple{typeof(block),EnergyForcesVirial}),
                 ("cell_task",PP.cell_task!,Tuple{typeof(w),typeof(p),EnergyForcesVirial,Int}))
        for (name,f,signature) in targets
            llvm=sprint(io -> code_llvm(io,f,signature;debuginfo=:none))
            native=sprint(io -> code_native(io,f,signature;debuginfo=:none,syntax=:intel))
            write(joinpath(directory,"$(name)-W$width.ll"),llvm)
            write(joinpath(directory,"$(name)-W$width.asm"),native)
            vecmul=count(line -> occursin(r"fmul <\d+ x double>",line),split(llvm,'\n'))
            packed=r"v(?:add|sub|mul|div|sqrt|fmadd|fmsub|fnmadd|fnmsub)[a-z0-9]*pd\s"
            zmmlines=filter(line -> occursin(packed,line) && occursin("zmm",line),split(native,'\n'))
            ymmlines=filter(line -> occursin(packed,line) && occursin("ymm",line),split(native,'\n'))
            push!(summary,"| $width | $name | $vecmul | $(length(zmmlines)) | $(length(ymmlines)) |")
            if !isempty(zmmlines)
                println("W=$width $name example: ",strip(first(zmmlines)))
            end
        end
    end
    append!(summary,["", "Counts are static code-inspection diagnostics, not executed instruction counts,",
        "a throughput measurement, or proof that the complete driver benefits from wider vectors.",
        "W=1 denotes one logical pair lane, not an absence of compiler auto-vectorisation:",
        "the compiler can emit vector instructions for other independent operations even at W=1.",
        "Callback arithmetic and fused traversal both need inspection: geometry packing, masks,",
        "lane occupancy, endpoint updates, reductions and scheduling also contribute to run time."])
    write(joinpath(directory,"README.md"),join(summary,'\n')*"\n")
    println(joinpath(directory,"README.md"))
end
inspect_widths()
