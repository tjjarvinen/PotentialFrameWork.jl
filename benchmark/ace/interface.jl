# julia --project=<application-env> --threads=4 benchmark/ace/interface.jl [output.toml]
using PotentialFrameWork, AtomsBase, Unitful, SpheriCart, LinearAlgebra, TOML, SIMD
const A=PotentialFrameWork.ACE
BLAS.set_num_threads(1)
model,ps,st=A.load_resolved(joinpath(@__DIR__,"../../data/ace_ref/resolved_model_v1.toml"))
sys=AtomsBase.atomic_system([:Ti=>[0.1,0.2,0.3]*u"Å",:Al=>[2.3,0.7,1.2]*u"Å",:Ti=>[1.1,2.4,1.5]*u"Å"],
    ([4.7,0.,0.]*u"Å",[0.5,4.8,0.]*u"Å",[0.1,0.3,4.6]*u"Å"),(true,true,true))
selection=A.BasisSelection(model;indices=1:20)
reduced=A.prune_basis(model,ps,selection)
cases=(("full",model,ps),("masked",A.mask_basis(model,selection),ps),
       ("pruned",reduced.model,reduced.parameters),
       ("species_0",A.project_species(model,[:Ti],0:0),ps),
       ("species_1_2",A.project_species(model,[:Ti],1:2),ps))
records=Dict{String,Any}()
for (label,m,p) in cases
    input=A.ace_inputs(sys,m;environment=A.FullEnvironment(width=4))
    env=input.environment; w=A.CPUExecutor(m,env;slots=4)
    execute=()->A.execute_prepared!(w,input.positions,input.cell,p)
    for _=1:4; execute(); end
    A.prepared_result(w); A.ace_inputs(sys,m;environment=env); A.prepare_executor!(w,env); execute()
    bytes=@allocated execute()
    times=map(1:7) do sample
        elapsed=@elapsed begin
            for _=1:100; execute(); end
        end
        elapsed/100
    end
    sort!(times)
    output_bytes=@allocated A.prepared_result(w)
    prepare_bytes=@allocated A.ace_inputs(sys,m;environment=env)
    A.prepare_executor!(w,env)
    records[label]=Dict("kernel_allocations_bytes"=>bytes,"seconds_min"=>first(times),
        "seconds_median"=>times[4],"seconds_max"=>last(times),
        "public_numeric_output_bytes"=>output_bytes,"boundary_preparation_bytes"=>prepare_bytes,
        "model_bytes"=>Base.summarysize(m),"parameters_bytes"=>Base.summarysize(p),
        "executor_bytes"=>Base.summarysize(w),"readout_cache_bytes"=>sizeof(w.readout),
        "precontracted_cache_bytes"=>sizeof(w.weights),"center_flags_bytes"=>sizeof(w.center_mask),
        "scratch_bytes_4_slots"=>Base.summarysize(w.scratch),
        "force_bytes"=>sizeof(w.force),"sites"=>length(sys),"active_edges"=>count(env.active),
        "invariants"=>size(m.coupling,1),"products"=>length(m.products),"density"=>length(m.density),
        "radial"=>length(m.radial_spec),"full_harmonics"=>(m.lmax+1)^2)
end
report=Dict("julia"=>string(VERSION),"threads"=>Threads.nthreads(),"blas_threads"=>BLAS.get_num_threads(),
    "cpu"=>Sys.cpu_info()[1].model,"SIMD"=>string(pkgversion(SIMD)),
    "SpheriCart"=>string(pkgversion(SpheriCart)),"Unitful"=>string(pkgversion(Unitful)),
    "method"=>"serial prepared E/F/W; 4 exclusive slots; 7 samples of 100 repetitions; full-system periodic three-atom environment; output/preparation allocations separate",
    "cases"=>records)
path=isempty(ARGS) ? joinpath(@__DIR__,"interface_measurements.toml") : first(ARGS)
open(io -> TOML.print(io,report;sorted=true),path,"w")
println("Wrote ",path)
