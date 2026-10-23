# Dedicated pinned producer; no fitting, and never writes the source recipe.
# julia --project=data/ace_ref data/ace_ref/export_resolved_model.jl
using ACEpotentials, SHA, TOML, LinearAlgebra, SparseArrays
const M = ACEpotentials.Models
VERSION == v"1.13.2" || error("Use Julia 1.13.2")
pkgversion(ACEpotentials) == v"0.10.2" || error("Unsupported producer version")

function plain_fields(obj)
    Dict(string(k)=>getproperty(obj,k) for k in propertynames(obj) if k != :meta)
end

function radial_export(basis)
    basis isa M.SplineRnlrzzBasis || error("Only resolved spline radial bases are supported")
    channels = map(CartesianIndices(basis.transforms)) do ij
        t=basis.transforms[ij]; env=basis.envelopes[ij]; spl=basis.splines[ij]
        t.trans isa M.GeneralizedAgnesiTransform || error("Unsupported radial transform")
        (env isa M.PolyEnvelope2sX || env isa M.ACE1_PolyEnvelope1sR) || error("Unsupported radial envelope")
        # The reconstructed object is a cubic spline, not the analytic polynomial
        # from which it was sampled. Export its exact piecewise cubic space.
        # Four evaluations determine each cubic; independent interior samples
        # below bound the floating-point change of coordinates.
        nodes=collect(spl.itp.ranges[1]); n=length(basis)
        coeff=zeros(4,n,length(nodes)-1)
        ts=[0.,1/3,2/3,1.]
        V=[t^q for t in ts,q=0:3]
        for k=1:length(nodes)-1
            lo=nodes[k]; step=nodes[k+1]-lo
            samples=permutedims(hcat([spl(lo+step*t) for t in ts]...))
            coeff[:,:,k] .= V \ samples
            for t in (0.137,0.513,0.927)
                out=vec(transpose([1.,t,t^2,t^3])*coeff[:,:,k])
                isapprox(out,spl(lo+step*t);atol=3e-13,rtol=3e-13) || error("Spline export failed")
            end
        end
        Dict("centre"=>ij[1],"neighbour"=>ij[2],
             "transform"=>merge(plain_fields(t.trans),Dict("yin"=>t.yin,"ycut"=>t.ycut,
                 "rcut"=>t.rcut,"normalization_rin"=>t.rin)),
             "envelope_kind"=>string(nameof(typeof(env))),"envelope"=>plain_fields(env),
             "nodes"=>nodes,"coefficient_shape"=>collect(size(coeff)),"coefficients"=>vec(coeff))
    end
    Dict("kind"=>"piecewise_cubic_agnesi", "spec"=>[Dict(string(k)=>v for (k,v) in pairs(b)) for b in basis.spec],
         "channels"=>vec(channels))
end

path=joinpath(@__DIR__,"reference_model.json")
hash=bytes2hex(sha256(read(path)))
p,saved=ACEpotentials.load_model(path); model=p.model
isempty(p.ps.rbasis) && isempty(p.ps.pairbasis) || error("Unexpected trainable radial parameters")
flat=vcat(vec(p.ps.WB),vec(p.ps.Wpair))
flat==saved["model_parameters"] || error("Raw parameter ordering mismatch")
length(flat)==270 || error("Unexpected parameter count")
model.tensor.LL==(0,) || error("Only scalar invariants supported")
model.tensor.aabasis.hasconst && error("Constant products need a new schema")
model.Vref isa M.OneBody || error("Only one-body fixed references supported")
C=model.tensor.A2Bmaps[1]; rows,cols,vals=findnz(C)
versions=Dict(name=>saved["manifest"]["deps"][name][1]["version"] for name in
    ("ACEpotentials","EquivariantTensors","SpheriCart","Polynomials4ML","Lux","ACEfit"))
resolved=Dict("schema"=>"PotentialFrameWork.resolved_ACE", "schema_version"=>1,
    "source_sha256"=>hash,"producer_julia"=>string(VERSION),"producer_versions"=>versions,
    "species_atomic_numbers"=>collect(model._i2z),
    "reference_eV"=>[model.Vref.E0[z] for z in model._i2z],
    "angular"=>Dict("kind"=>"real_spherical", "normalization"=>"L2", "ordering"=>"l*l+l+m+1",
        "lmax"=>isqrt(length(model.ybasis))-1),
    "radial"=>radial_export(model.rbasis),"pair_radial"=>radial_export(model.pairbasis),
    "density_indices"=>[collect(t) for t in model.tensor.abasis.spec],
    "products"=>[collect(t) for spec in model.tensor.aabasis.specs for t in spec],
    "coupling"=>Dict("shape"=>collect(size(C)),"rows"=>rows,"columns"=>cols,"values"=>vals),
    "parameters"=>Dict("WB_shape"=>collect(size(p.ps.WB)),"Wpair_shape"=>collect(size(p.ps.Wpair)),
        "flat"=>flat,"ordering"=>"vec(WB), vec(Wpair); Julia column-major"))
open(joinpath(@__DIR__,"resolved_model_v1.toml"),"w") do io
    TOML.print(io,resolved;sorted=true)
end
bytes2hex(sha256(read(path)))==hash || error("Source model changed")
println("Exported 270 parameters, ",length(model.tensor.abasis)," densities, ",size(C)," coupling map")
