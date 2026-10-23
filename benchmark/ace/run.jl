# Run without concurrent tests/benchmarks:
# julia --project=data/ace_ref --threads=4 benchmark/ace/run.jl
using ACEpotentials
pushfirst!(LOAD_PATH,normpath(joinpath(@__DIR__,"../..")))
using PotentialFrameWork, ACEpotentials, ExtXYZ, AtomsCalculators, Unitful
using LinearAlgebra, TOML
using StaticArrays
using AtomsBase
const A=PotentialFrameWork.ACE
BLAS.set_num_threads(1)

function measurement(f;samples=9,warmup=3)
    for _=1:warmup; f(); end
    seconds=Float64[]; bytes=Int[]
    for _=1:samples
        r=@timed f(); push!(seconds,r.time); push!(bytes,r.bytes)
    end
    sorted=sort(seconds)
    Dict("seconds"=>seconds,"bytes"=>bytes,"minimum_s"=>first(sorted),
         "median_s"=>sorted[cld(samples,2)],"maximum_s"=>last(sorted),"median_bytes"=>sort(bytes)[cld(samples,2)])
end
function generated_angular!(Y,dY,basis,Rs,gradients)
    # A function barrier preserves the concrete generated-basis type across the
    # neighbour loop; dispatching a runtime-L basis once per neighbour would
    # measure boxing/dispatch rather than its generated numerical kernel.
    for i in eachindex(Rs)
        if gradients
            y,dy=ACEpotentials.Models.P4ML.SpheriCart.compute_with_gradients(basis,Rs[i])
            for n in eachindex(y); Y[i,n]=y[n]; dY[i,n]=dy[n]; end
        else
            y=ACEpotentials.Models.P4ML.SpheriCart.compute(basis,Rs[i])
            for n in eachindex(y); Y[i,n]=y[n]; end
        end
    end
    nothing
end
function run_benchmarks()
refdir=normpath(joinpath(@__DIR__,"../../data/ace_ref"))
model,ps,st=A.load_recipe(joinpath(refdir,"reference_model.json"))
oracle,_=ACEpotentials.load_model(joinpath(refdir,"reference_model.json"))
frames=ExtXYZ.load(joinpath(refdir,"reference_results.xyz"))
systems=Tuple{Int,Any}[(1,frames[1]),(33,frames[33])]
# Additional controlled neighbour-count and species-group workloads. Timings are
# specific to the supplied basis; only the geometries/compositions vary here.
for g in (2,4),spacing in (2.8,4.0),species_count in (1,2)
    positions=vec([SVector(x,y,z)*spacing for x=0:g-1,y=0:g-1,z=0:g-1])
    labels=[species_count==1 || isodd(i) ? :Ti : :Al for i in eachindex(positions)]
    H=ntuple(a -> SVector(ntuple(b -> a==b ? g*spacing : 0.,3))*u"Å",3)
    sys=AtomsBase.atomic_system([labels[i]=>positions[i]*u"Å" for i in eachindex(labels)],H,(true,true,true))
    push!(systems,(-length(systems),sys))
end
records=Dict[]
for (frame,sys) in systems
    for width in (2,4,8),backend in (:sphericart,:generated,:block),fused in (false,true)
        env=A.FullEnvironment(;width); input=A.ace_inputs(sys,model;environment=env)
        w=A.CPUExecutor(model,env;slots=Threads.nthreads(),angular_backend=backend)
        ref=A.ace_evaluate(input.positions,model,ps,env,input.cell)
        A.execute_prepared!(w,input.positions,input.cell,ps;parallel=true,fused)
        checked=A.prepared_result(w)
        isapprox(checked.forces,ref.forces;atol=1e-9,rtol=1e-10) || error("benchmark force mismatch")
        prepared=measurement(() -> A.execute_prepared!(w,input.positions,input.cell,ps;parallel=true,fused))
        complete=measurement() do
            inp=A.ace_inputs(sys,model;environment=env)
            A.prepare_executor!(w,env)
            A.execute_prepared!(w,inp.positions,inp.cell,ps;parallel=true,fused)
            A.prepared_result(w)
        end
        push!(records,Dict("frame"=>frame,"atoms"=>length(sys),"width"=>width,"angular"=>string(backend),"fused"=>fused,
            "prepared_EFW"=>prepared,"complete_EFW"=>complete,"force_bytes"=>sizeof(w.force),
            "scratch_bytes"=>Base.summarysize(w.scratch),"model_bytes"=>Base.summarysize(model),
            "edges"=>count(env.active),"packed_lanes"=>length(env.active),"phase_sizes"=>length.(env.phases),
            "species_counts"=>[count(==(z),env.species) for z in eachindex(model.atomic_numbers)],
            "preparation_stage_seconds_last_sample"=>Dict(string(k)=>v for (k,v) in env.timings)))
    end
    env=A.FullEnvironment(width=4); inp=A.ace_inputs(sys,model;environment=env)
    calc=A.ACECalculator(model,ps;slots=Threads.nthreads())
    data=A.site_diagnostics(inp.positions,model,ps,env,inp.cell,1)
    species=[env.species[env.neighbours[k]] for k in data.edges]
    angular_basis=A.angular_basis(model.lmax,Float64)
    push!(records,Dict("frame"=>frame,"atoms"=>length(sys),"boundary_measurements"=>Dict(
        "preparation"=>measurement(() -> A.ace_inputs(sys,model;environment=env)),
        "public_EFW"=>measurement(() -> AtomsCalculators.energy_forces_virial(sys,calc;parallel=true)),
        "public_energy"=>measurement(() -> AtomsCalculators.potential_energy(sys,calc)),
        "serial_independent_EFW"=>measurement(() -> A.ace_evaluate(inp.positions,model,ps,env,inp.cell)),
        "oracle_public_EFW"=>measurement(() -> AtomsCalculators.energy_forces_virial(sys,oracle)),
        "first_site_radial_values_derivatives"=>measurement(() -> [A.radial_ed(model.radial[env.species[1],species[j]],data.rs[j]) for j in eachindex(species)]),
        "first_site_angular_values_derivatives"=>measurement(() -> A.angular_ed(angular_basis,data.Rs)),
        "first_site_reverse"=>measurement(() -> A.site_reverse(model,ps,env.species[1],data)),
        "design_EFW"=>measurement(() -> A.design_matrix(inp.positions,model,ps,env,inp.cell)))))
    if frame==1
        Z=ACEpotentials.Models.Zygote
        params=A.trainable_parameters(model,ps)
        loss=(x,p) -> sum(sin.(inp.positions).*A.ace_evaluate(x,model,p,env,inp.cell).forces)
        _,back=Z.pullback(loss,inp.positions,params)
        push!(records,Dict("frame"=>frame,"AD_force_loss"=>Dict(
            "forward_and_save"=>measurement(() -> Z.pullback(loss,inp.positions,params)),
            "saved_backward"=>measurement(() -> back(1.)),
            "saved_pullback_bytes"=>Base.summarysize(back))))
    end
end

# Standalone angular workloads. Explicit blocks and generated scalar functions
# are measured as implemented, including their actual temporary allocations.
# Block value-only still uses its joint value/gradient recurrence internally.
angular=Dict[]
SC=ACEpotentials.Models.P4ML.SpheriCart
for L in (1,2,4,6),n in (16,64,256),width in (2,4,8),gradients in (false,true)
    Rs=[SVector(sin(i+0.3),cos(0.7i),sin(2i+0.9)) for i=1:n]
    Y=zeros(n,(L+1)^2); dY=zeros(SVector{3,Float64},size(Y))
    block=A.BlockAngular(L;width)
    batched=SC.SphericalHarmonics(L;static=false)
    generated=SC.SphericalHarmonics(L;static=true)
    block_time=measurement(() -> A.angular_block!(Y,gradients ? dY : nothing,block,Rs))
    batch_time=measurement(() -> gradients ? SC.compute_with_gradients!(Y,dY,batched,Rs) : SC.compute!(Y,batched,Rs))
    generated_time=measurement(() -> generated_angular!(Y,dY,generated,Rs,gradients))
    push!(angular,Dict("lmax"=>L,"neighbours"=>n,"width"=>width,"gradients"=>gradients,
        "sphericart_batched"=>batch_time,"sphericart_generated_single"=>generated_time,"explicit_blocks"=>block_time))
end
report=Dict("julia"=>string(VERSION),"threads"=>Threads.nthreads(),"blas_threads"=>BLAS.get_num_threads(),
    "cpu"=>Sys.cpu_info()[1].model,"cpu_target"=>get(ENV,"JULIA_CPU_TARGET","native/default"),
    "affinity"=>"not pinned", "model_sha256"=>model.source_sha256,"samples"=>9,"warmup"=>3,
    "versions"=>Dict(pkg.name=>string(pkgversion(mod)) for (pkg,mod) in Base.loaded_modules
        if pkg.name in ("ACEpotentials","ACEfit","SpheriCart","Lux","LuxCore","Zygote","ForwardDiff","SIMD","StaticArrays","AtomsBase","AtomsCalculators","EquivariantTensors","Polynomials4ML","Unitful")),
    "full_stack"=>records,"angular"=>angular)
path=joinpath(@__DIR__,"results-julia$(VERSION)-t$(Threads.nthreads()).toml")
open(path,"w") do io; TOML.print(io,report;sorted=true); end
println(path)
end
run_benchmarks()
