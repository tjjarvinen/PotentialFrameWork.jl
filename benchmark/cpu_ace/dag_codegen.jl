# Inspect concrete production callees, not only the dispatch/wrapper body.
include("common.jl")
function main()
    out=abspath(ARGS[1]); mkpath(out); records=Dict[]
    for T in (Float32,Float64)
        m,ps,_=A.load_resolved(joinpath(@__DIR__,"results/2026-11-10/models/large-tial-resolved.toml");precision=T)
        sys=structure(bcc_geometry(4,m.atomic_numbers);T)
        input=A.ace_inputs(sys,m;environment=A.FullEnvironment(width=4),precision=T)
        w=A.CPUExecutor(m,input.environment;precision=T,product_algorithm=:dag)
        A.execute_prepared!(w,input.positions,input.cell,ps)
        s=w.scratch[1]; plan=w.prepared.product_plan
        site=NL.SiteView(input.environment,input.positions,input.cell,1)
        A.NI.evaluate_site!(w.slots[1].outputs[1],w.prepared,site,ps,s,(AC.Energy(),AC.Forces(),AC.Virial()))
        weights=view(w.weights,:,A.NI.central_species(site))
        jobs=(("forward",A.dag_forward!,(s.P,s.U,s.V,plan,s.dag_values)),
              ("reverse",A.dag_reverse!,(s.dA,s.dV,weights,plan,s.dag_values,s.dag_adjoints)),
              ("combined",A.dag_products!,(s.P,s.dA,s.dV,s.U,s.V,weights,plan,s.dag_values,s.dag_adjoints)))
        for (name,f,args) in jobs
            f(args...)
            bytes=@allocated f(args...)
            @assert bytes==0
            prefix=joinpath(out,name*"-"*string(T))
            open(prefix*".llvm.txt","w") do io; code_llvm(io,f,typeof(args)); end
            open(prefix*".native.txt","w") do io; code_native(io,f,typeof(args)); end
            open(prefix*".warntype.txt","w") do io; code_warntype(io,f,typeof(args)); end
            push!(records,Dict("callee"=>name,"precision"=>string(T),"host_bytes"=>bytes,"signature"=>string(typeof(args))))
        end
    end
    open(io->TOML.print(io,Dict("julia"=>string(VERSION),"script_sha256"=>bytes2hex(sha256(read(@__FILE__))),"common_sha256"=>bytes2hex(sha256(read(joinpath(@__DIR__,"common.jl")))),"source_files_sha256"=>source_provenance(dirname(dirname(pathof(PotentialFrameWork)))),"callees"=>records);sorted=true),joinpath(out,"inspection.toml"),"w")
end
main()
