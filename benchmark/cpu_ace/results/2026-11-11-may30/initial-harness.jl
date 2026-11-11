# Dedicated live integration/benchmark. See MAY30-REPORT.md and may30_setup.jl.
using ACEpotentials, PotentialFrameWork, SpheriCart, BenchmarkTools
using AtomsBase, AtomsCalculators, Unitful, StaticArrays
using LinearAlgebra, SparseArrays, SHA, TOML, Statistics, Pkg
include(joinpath(@__DIR__,"../../data/ref_model/train_may30_ace.jl"))
include(joinpath(@__DIR__,"../../data/ace_ref/resolved_export.jl"))
include(joinpath(@__DIR__,"may30_fixtures.jl"))
const A=PotentialFrameWork.ACE
const AC=AtomsCalculators
BLAS.set_num_threads(1)
const OUT=abspath(get(ARGS,1,joinpath(@__DIR__,"results/2026-11-11-may30")))
mkpath(OUT)
write_toml(path,d)=open(io->TOML.print(io,d;sorted=true),path,"w")
hashfile(p)=bytes2hex(sha256(read(p)))
function system(d)
    atomic_system([ChemicalSpecies(d.species[i])=>SVector{3}(d.positions[:,i])*u"Å" for i in axes(d.positions,2)],
        ntuple(j->SVector{3}(d.cell[:,j])*u"Å",3),(true,true,true))
end
function numeric(y)
    E=ustrip(u"eV",y.energy); F=reduce(hcat,[ustrip.(u"eV/Å",v) for v in y.forces]); W=Matrix(ustrip.(u"eV",y.virial))
    isfinite(E) && all(isfinite,F) && all(isfinite,W) || error("nonfinite E/F/V")
    size(W)==(3,3) && size(F,1)==3 || error("shape mismatch")
    (;energy=E,forces=F,virial=W)
end
function compare(x,y,n;approximation=false)
    size(x.forces)==size(y.forces)==(3,n) && size(x.virial)==size(y.virial)==(3,3) || error("shape mismatch")
    errors=[abs(x.energy-y.energy)/n,maximum(abs,x.forces-y.forces),maximum(abs,x.virial-y.virial)]
    # Parity uses per-atom energy without Eref-scaled relative slack.
    absolute=approximation ? [2e-6,2e-5,2e-4] : [2e-10,2e-9,2e-8]
    scales=[1.,maximum(abs,x.forces),maximum(abs,x.virial)]
    rel=approximation ? 2e-7 : 2e-11
    thresholds=absolute+rel*scales
    all(errors .<= thresholds) || error("parity failed: $errors > $thresholds (approximation=$approximation)")
    Dict("energy_per_atom_eV"=>errors[1],"force_component_eV_per_A"=>errors[2],"virial_component_eV"=>errors[3],"thresholds"=>thresholds,"passed"=>true)
end
function measure(f)
    warm=@elapsed f()
    samples=parse(Int,get(ENV,"PFW_MAY30_SAMPLES","10")); seconds=parse(Float64,get(ENV,"PFW_MAY30_SECONDS","10"))
    t=run(@benchmarkable($f(),evals=1);samples,seconds)
    Dict("raw_ns"=>t.times,"gc_ns"=>t.gctimes,"median_ns"=>median(t).time,"min_ns"=>minimum(t.times),"max_ns"=>maximum(t.times),
        "bytes"=>t.memory,"allocations"=>t.allocs,"samples"=>length(t.times),"undersampled"=>length(t.times)<10,"warmup_s"=>warm)
end
checkpoint=joinpath(@__DIR__,"../../data/ref_model/may30_d10_n500.jls")
construction=Dict{String,Any}()
construction["reload_s"]=@elapsed calc,saved=load_may30_fit(checkpoint)
@assert saved.settings.degree==10 && saved.settings.order==3 && saved.settings.rcut==6 && length(saved.train_indices)==500
@assert collect(calc.model._i2z)==[22,13,23]
construction["et_conversion_s"]=@elapsed et=ACEpotentials.ETModels.convert2et_full(calc.model,calc.ps,calc.st;rng=MersenneTwister(1234))
@assert length(et.calcs)==3 && [c.rcut for c in et.calcs]==[0.,6.,6.]
et_species=[[Int(atomic_number(z)) for z in et.calcs[2].st.rembed.rbasis.trans.zlist],
            [Int(atomic_number(z)) for z in et.calcs[3].st.rembed.trans.zlist]]
@assert all(==([22,13,23]),et_species)
nnodes=parse(Int,get(ENV,"PFW_MAY30_NODES","1000"))
construction["splinify_s"]=@elapsed sm=M.splinify(calc.model,calc.ps;nnodes)
# The v0.10.3 calculator-level splinify loses the fitted readouts. Use model API.
_,sst=Lux.setup(MersenneTwister(1234),sm)
sps=merge(calc.ps,(rbasis=NamedTuple(),pairbasis=NamedTuple()))
spl=ACEpotentials.ACEPotential(sm,sps,sst)
@assert spl.ps.WB==calc.ps.WB && spl.ps.Wpair==calc.ps.Wpair && sm.Vref===calc.model.Vref
resolved=calculator_resolved(spl;source_sha256=hashfile(checkpoint))
artifact=joinpath(OUT,"may30-resolved.toml"); write_toml(artifact,resolved)
construction["pfw_import_s"]=@elapsed m,ps,st=A.load_resolved(artifact)
A.ProductDAG(m.products,length(m.density)) # warm plan compilation outside setup measurement
plan_trial=@timed A.ProductDAG(m.products,length(m.density))
construction["warmed_dag_plan_s"]=plan_trial.time
construction["warmed_dag_plan_allocated_bytes"]=plan_trial.bytes
versions=Dict(info.name=>string(info.version) for info in values(Pkg.dependencies()) if !isnothing(info.version))
meta=Dict("starting_commit"=>"fbbe8c4efb96ad935d2c006b0b1971998b832748","julia"=>string(VERSION),"versions"=>versions,
    "checkpoint_sha256"=>hashfile(checkpoint),"training_script_sha256"=>hashfile(joinpath(dirname(checkpoint),"train_may30_ace.jl")),
    "checkpoint_julia"=>string(saved.julia_version),"checkpoint_ACEpotentials"=>string(saved.acepotentials_version),
    "settings"=>Dict(string(k)=>(v isa Tuple ? collect(v) : v) for (k,v) in pairs(saved.settings)),"training_frames"=>length(saved.train_indices),"solver"=>saved.solver,"smoothness"=>saved.smoothness,
    "fit_weights"=>saved.weights,"train_indices"=>saved.train_indices,"heldout_indices"=>saved.heldout_indices,
    "source_data_sha256"=>saved.data_sha256,"resolved_sha256"=>hashfile(artifact),"spline_nodes"=>nnodes,
    "types"=>Dict("polynomial"=>string(typeof(calc)),"et"=>string(typeof(et)),"spline"=>string(typeof(spl))),
    "angular"=>resolved["angular"],"et_cutoffs_A"=>[c.rcut for c in et.calcs],"et_species"=>et_species,
    "width"=>4,"simd"=>true,"parallel"=>false,"fusion"=>false,"angular_backend"=>"sphericart","species"=>m.atomic_numbers,"cutoff_A"=>6.,"WB_shape"=>collect(size(ps.WB)),"Wpair_shape"=>collect(size(ps.Wpair)),
    "density_count"=>length(m.density),"dag_multiply_nodes"=>length(plan_trial.value.left),
    "product_orders"=>Dict(string(k)=>count(p->length(p)==k,m.products) for k in sort(unique(length.(m.products)))),"products"=>length(m.products),"coupling_shape"=>collect(size(m.coupling)),"coupling_nnz"=>nnz(m.coupling),
    "radial_count"=>length(m.radial_spec),"pair_radial_count"=>length(m.pair_spec),"threads"=>Threads.nthreads(),"blas_threads"=>BLAS.get_num_threads(),
    "affinity"=>read(`/usr/bin/taskset -pc $(getpid())`,String),"construction"=>construction)
write_toml(joinpath(OUT,"metadata.toml"),meta)
provenance=Dict{String,Any}()
for info in values(Pkg.dependencies())
    info.name in ("ACEpotentials","EquivariantTensors","Polynomials4ML","SpheriCart","PotentialFrameWork") || continue
    files=Dict{String,String}()
    for (dir,_,names) in walkdir(joinpath(info.source,"src")),name in names
        file=joinpath(dir,name); files[relpath(file,info.source)]=hashfile(file)
    end
    provenance[info.name]=Dict("version"=>string(info.version),"tree_hash"=>string(info.tree_hash),"source"=>info.source,"files_sha256"=>files)
end
write_toml(joinpath(OUT,"source-provenance.toml"),provenance)
# Independent angular and radial source samples, including radial dependence.
Rs=[SVector{3}(r.*v/norm(v)) for r in (.7,1.9,4.1),v in ([.3,.7,-.2],[-.5,.1,.9])][:]
sourceY,sourceDY=SpheriCart.compute_with_gradients(calc.model.ybasis.scbasis,Rs)
for backend in (:sphericart,:generated,:block)
    b=A.make_angular(m,Float64,backend,4); Y=similar(sourceY); dY=similar(sourceDY)
    A.angular_ed!(Y,dY,b,Rs)
    @assert isapprox(Y,sourceY;atol=2e-11,rtol=2e-12) && isapprox(dY,sourceDY;atol=2e-11,rtol=2e-12)
end
radial_records=Dict[]
for (name,source,channels) in (("radial",sm.rbasis,m.radial),("pair",sm.pairbasis,m.pair_radial)),iz=1:3,jz=1:3
    ev=0.; ed=0.
    for r in (1.,2.3,3.1,5.9,6.0 - 1e-6,6.,6.0 + 1e-6)
        ch=channels[iz,jz]; x=A.transformed(ch.transform,r)
        env=source.envelopes[iz,jz]
        @assert A.envelope(ch,r,x) ≈ M.evaluate(env,r,x)
        value,deriv=A.radial_ed(ch,r)
        # Differentiate source's actual spline, transform and envelope jointly.
        f(t)=source.splines[iz,jz](M.evaluate(source.transforms[iz,jz],t))*M.evaluate(env,t,M.evaluate(source.transforms[iz,jz],t))
        sv=f(r); sd=r==6.0 ? zero(sv) : ACEpotentials.Models.ForwardDiff.derivative(f,r)
        # At exactly rcut both envelopes vanish with zero slope. The source
        # Interpolations Dual bounds check rejects the endpoint with nonzero
        # transform tangent, so use this analytic endpoint limit only there.
        ev=max(ev,maximum(abs,value-sv)); ed=max(ed,maximum(abs,deriv-sd))
    end
    @assert ev<2e-9 && ed<2e-8
    push!(radial_records,Dict("family"=>name,"centre"=>iz,"neighbour"=>jz,"value_error"=>ev,"derivative_error"=>ed))
end
write_toml(joinpath(OUT,"basis-validation.toml"),Dict("angular_passed"=>true,"radial"=>radial_records))
records=Dict[]
sizes=parse.(Int,split(get(ENV,"PFW_MAY30_SIZES","4,8,12"),','))
fixtures=vcat([may30_geometry(:bcc;g) for g in sizes],[may30_geometry(:hcp),may30_geometry(:omega)])
for (index,d) in enumerate(fixtures)
    name=d.kind*"-"*string(length(d.species)); println("Validating ",name); flush(stdout)
    sys=system(d); n=length(sys); H=d.cell; X=d.positions
    heights=[1/norm(inv(H)[a,:]) for a=1:3]; @assert minimum(heights)>12
    @assert Set(d.species)==Set([13,22,23]) && d.species[1]==d.species[2]==23
    @assert 0<norm(X[:,2]-X[:,1])<6
    inp=A.ace_inputs(sys,m;environment=A.FullEnvironment(width=4)); env=inp.environment
    distances=[norm(A.edge_displacement(X,H,env,i,k)) for i=1:n for k=env.offsets[1,i]:env.offsets[end,i]-1 if env.active[k]]
    @assert minimum(distances)>1.5
    nearest1=minimum(norm(A.edge_displacement(X,H,env,1,k)) for k=env.offsets[1,1]:env.offsets[end,1]-1 if env.active[k])
    @assert isapprox(nearest1,d.vv_distance;atol=1e-12,rtol=1e-12)
    @assert d.species[d.vv_neighbor]==23 && d.vv_distance<6
    jobs=["polynomial"=>calc,"et_polynomial"=>et,"spline"=>spl]
    ys=Dict(key=>numeric(AC.energy_forces_virial(sys,c)) for (key,c) in jobs)
    @assert maximum(abs,ys["polynomial"].forces)>1e-5 && maximum(abs,ys["polynomial"].virial)>1e-4
    checks=Dict("ET_conversion"=>compare(ys["polynomial"],ys["et_polynomial"],n),
        "spline_approximation"=>compare(ys["polynomial"],ys["spline"],n;approximation=true))
    executors=Dict{String,Any}(); public=Dict{String,Any}()
    setup=Dict{String,Any}()
    for algorithm in (:optimized,:dag)
        key=String(algorithm)
        setup[key*"_calculator_s"]=@elapsed c=A.ACECalculator(m,ps,st;slots=1,product_algorithm=algorithm)
        setup[key*"_executor_s"]=@elapsed ex=A.CPUExecutor(m,env;slots=1,product_algorithm=algorithm)
        A.execute_prepared!(ex,X,H,ps;parallel=false,simd=true,fused=false,product_algorithm=algorithm)
        checks[key*"_prepared"]=compare(ys["spline"],A.prepared_result(ex),n)
        checks[key*"_public"]=compare(ys["spline"],numeric(AC.energy_forces_virial(sys,c)),n)
        executors[key]=ex; public[key]=c
        @assert A.effective_product_algorithm(algorithm,m.policy.orders)==algorithm
    end
    # Save full precision exact geometries/oracles for runtime-independent regression.
    fixture=Dict("name"=>name,"positions"=>vec(X),"cell"=>vec(H),"species"=>d.species,
        "energy"=>ys["spline"].energy,"forces"=>vec(ys["spline"].forces),"virial"=>vec(ys["spline"].virial))
    write_toml(joinpath(OUT,name*"-fixture.toml"),fixture)
    timings=Dict{String,Any}()
    if get(ENV,"PFW_MAY30_TIMING","true")=="true"
        for (key,c) in jobs; timings[key*"_public"]=measure(()->AC.energy_forces_virial(sys,c)); end
        for algorithm in (:optimized,:dag)
            key=String(algorithm); ex=executors[key]; c=public[key]
            timings[key*"_prepared"]=measure(()->A.execute_prepared!(ex,X,H,ps;parallel=false,simd=true,fused=false,product_algorithm=algorithm))
            timings[key*"_public"]=measure(()->AC.energy_forces_virial(sys,c;parallel=false,simd=true,fused=false,product_algorithm=algorithm))
        end
    end
    push!(records,Dict("fixture"=>name,"atoms"=>n,"fixture_sha256"=>hashfile(joinpath(OUT,name*"-fixture.toml")),
        "geometry_sha256"=>bytes2hex(sha256(string(X,H,d.species))),
        "seed"=>d.seed,"cell"=>vec(H),"composition"=>[count(==(z),d.species) for z in m.atomic_numbers],
        "face_heights_A"=>heights,"minimum_distance_A"=>minimum(distances),"VV_distance_A"=>d.vv_distance,"VV_neighbor"=>d.vv_neighbor,"VV_image"=>collect(d.vv_shift),
        "directed_neighbors"=>length(distances),"checks"=>checks,"timings"=>timings,"setup"=>setup))
    write_toml(joinpath(OUT,"records.toml"),Dict("records"=>records))
    println(name," passed; ",[(k,round(v["median_ns"]/1e6;digits=3)) for (k,v) in timings]); flush(stdout)
end
