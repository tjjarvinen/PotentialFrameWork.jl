using PotentialFrameWork, SpheriCart, BenchmarkTools, TOML, SHA, LinearAlgebra, Statistics, InteractiveUtils, SparseArrays
const A=PotentialFrameWork.ACE
const NL=PotentialFrameWork.NeighbourList
const AB=A.AtomsBase
const U=A.Unitful
const AC=A.AC
include(joinpath(@__DIR__,"../gpu_ace/workloads.jl"))
include(joinpath(@__DIR__,"../gpu_ace/provenance.jl"))
function structure(d;T=Float64)
    AB.atomic_system([AB.ChemicalSpecies(d.species[i])=>A.SVector{3}(T.(d.positions[:,i]))*U.u"Å" for i=1:size(d.positions,2)],
        ntuple(a->A.SVector{3}(T.(d.cell[:,a]))*U.u"Å",3),(true,true,true))
end
function measure(f,samples,seconds)
    first=@elapsed f()
    t=run(@benchmarkable($f(),evals=1);samples,seconds)
    times=t.times*1e-9
    Dict("median_s"=>median(times),"min_s"=>minimum(times),"max_s"=>maximum(times),"raw_s"=>times,
        "raw_gc_ns"=>t.gctimes,"samples"=>length(times),"undersampled"=>length(times)<10,
        "host_bytes"=>t.memory,"host_allocations"=>t.allocs,"warmup_s"=>first)
end
function accuracy(y,r)
    Dict("energy_abs_eV"=>abs(Float64(y.energy)-r.energy),
        "forces_max_eV_per_A"=>maximum(abs,Float64.(y.forces)-r.forces;init=0.),
        "virial_max_eV"=>maximum(abs,Float64.(y.virial)-r.virial;init=0.))
end
function agree(y,r,T)
    tol=T===Float32 ? 3e-4 : 2e-8
    @assert isapprox(y.energy,r.energy;atol=tol,rtol=tol)
    @assert isapprox(y.forces,r.forces;atol=tol,rtol=tol)
    @assert isapprox(y.virial,r.virial;atol=tol,rtol=tol)
end
function metadata(m,ps)
    orders=Dict(string(k)=>count(p->length(p)==k,m.products) for k in sort(unique(length.(m.products))))
    Dict("species"=>m.atomic_numbers,"radial_count"=>length(m.radial_spec),"pair_radial_count"=>length(m.pair_spec),
        "harmonic_count"=>(m.lmax+1)^2,"lmax"=>m.lmax,"density_count"=>length(m.density),"products"=>length(m.products),
        "product_orders"=>orders,"coupling_shape"=>collect(size(m.coupling)),"coupling_nnz"=>nnz(m.coupling),
        "active_readout_per_central_species"=>[count(!iszero,ps.WB[:,z]) for z in axes(ps.WB,2)],
        "active_product_weights_per_central_species"=>[count(!iszero,transpose(m.coupling)*ps.WB[:,z]) for z in axes(ps.WB,2)],
        "cutoff_A"=>A.NI.search_cutoff(m),"radial_representation"=>string(typeof(first(m.radial))),
        "coordinate_kind"=>string(m.coordinate_kind))
end
# Isolated microbenchmarks operate on warmed first-site intermediates. They
# deliberately do not introduce timers into production loops, nor sum to a total.
function diagnostics(w,X,H,ps,algorithm,samples,seconds,codepath)
    A.execute_prepared!(w,X,H,ps;product_algorithm=algorithm)
    m=w.model; s=w.scratch[1]; env=w.environment; site=NL.SiteView(env,X,H,1)
    n=A.NI.neighbour_count(site); iz=A.NI.central_species(site)
    props=(AC.Energy(),AC.Forces(),AC.Virial()); out=w.slots[1].outputs[1]
    p=typeof(w.prepared)(m,w.basis,w.weights,true,false,algorithm)
    A.NI.evaluate_site!(out,p,site,ps,s,props)
    # Capture first-site rather than last-site density state.
    productargs=(s.P,s.dA,s.dV,s.U,s.V,m.products,view(w.weights,:,iz),m.policy.orders,s.polynomial,s.reverse_polynomial,algorithm)
    jobs=Pair{String,Any}[
        "geometry_radii"=>(() -> [norm(A.NI.displacement(site,j)) for j=1:n]),
        "angular"=>(() -> A.angular_ed!(view(s.Y,1:n,:),view(s.dY,1:n,:),w.basis,view(s.Rs,1:n))),
        "density"=>(function ()
            fill!(s.A,0); fill!(s.U,0); fill!(s.V,0)
            for z in eachindex(m.atomic_numbers)
                r=A.NI.species_range(site,z); A.project_rows!(s,m,first(r),last(r),z in m.policy.species)
            end
        end),
        "products_reverse"=>(() -> A.prepared_products!(productargs...)),
        "readout"=>(() -> dot(view(w.weights,:,iz),s.P)),
        "site_math"=>(() -> A.NI.evaluate_site!(out,p,site,ps,s,props)),
        "reduction"=>(function ()
            NL.reduce_site!(w.force,w.inverse,w.energies,w.virials,env,1,site,out,props)
        end),
        "owning_output"=>(() -> A.prepared_result(w))]
    push!(jobs,"radial"=>(function ()
        for z in eachindex(m.atomic_numbers)
            r=A.NI.species_range(site,z)
            for k=first(r):env.width:last(r)
                rs=ntuple(lane->norm(s.Rs[k+lane-1]),env.width)
                active=ntuple(lane->A.NI.active(site,k+lane-1),env.width)
                A.radial_block!(s.R,s.dR,rs,m.radial[iz,z],k,Val(env.width),active,A.radial_parameters(ps,:custom_radial))
                A.radial_block!(s.Q,s.dQ,rs,m.pair_radial[iz,z],k,Val(env.width),active,A.radial_parameters(ps,:custom_pair))
            end
        end
    end))
    # Exact existing edge back-contraction, without reduction, at first-site state.
    push!(jobs,"force_back_contraction"=>(function ()
        for j=1:n
            A.NI.active(site,j) || continue
            d=s.Rs[j]; direction=d/norm(d); g=zero(d)
            for (a,(r,h)) in enumerate(m.density)
                adj=A.NI.neighbour_species(site,j) in m.policy.species ? s.dV : s.dA
                g+=adj[a]*(s.dR[j,r]*s.Y[j,h]*direction+s.R[j,r]*s.dY[j,h])
            end
            for q=1:size(s.Q,2); g+=ps.Wpair[q,iz]*s.dQ[j,q]*direction; end
            A.NI.set_gradient!(out,j,g)
        end
    end))
    results=Dict{String,Any}()
    for (name,f) in jobs
        A.NI.evaluate_site!(out,p,site,ps,s,props)
        results[name]=measure(f,samples,seconds)
    end
    if !isempty(codepath)
        mkpath(codepath)
        ch=m.radial[iz,1]; W=env.width
        args=(s.R,s.dR,ntuple(_->eltype(X)(2.3),W),ch,1,Val(W),ntuple(_->true,W),A.radial_parameters(ps,:custom_radial))
        for (name,f,args) in (("radial",A.radial_block!,args),("products_dispatch",A.prepared_products!,productargs),
                ("products_unrestricted",A.unrestricted_products!,(productargs[1:7]...,s.polynomial)),
                ("products_reference",A.product_projection!,productargs[1:10]),("site",A.NI.evaluate_site!,(out,p,site,ps,s,props)))
            open(joinpath(codepath,name*".llvm.txt"),"w") do io; code_llvm(io,f,typeof(args)); end
            open(joinpath(codepath,name*".native.txt"),"w") do io; code_native(io,f,typeof(args)); end
        end
    end
    results
end
