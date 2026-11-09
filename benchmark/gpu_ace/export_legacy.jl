# Isolated producer only: julia --project=benchmark/gpu_ace/legacy export_legacy.jl SOURCE.json OUTPUT.toml
using ACEpotentials, JSON, TOML, SHA, LinearAlgebra, SparseArrays
const J=ACEpotentials.JuLIP
const Old=ACEpotentials.ACE1
include(joinpath(@__DIR__,"../../src/submodules/ACE/legacy_radial.jl"))
pkgversion(ACEpotentials)==v"0.6.6" || error("Use isolated ACEpotentials 0.6.6")
encode(x::NamedTuple)=Dict("kind"=>"namedtuple","names"=>string.(collect(keys(x))),"items"=>[encode(v) for v in values(x)])
encode(x::AbstractArray)=Dict("kind"=>"array","eltype"=>string(eltype(x)),"shape"=>collect(size(x)),"values"=>collect(vec(x)))
encode(x::Number)=Dict("kind"=>"scalar","eltype"=>string(typeof(x)),"value"=>x)
function config(radial,z,z0,offset,cutoff)
    trans=radial.trans
    hasproperty(trans,:transforms) && (trans=trans.transforms[Old.z2i(trans.zlist,z),Old.z2i(trans.zlist,z0)])
    affine=hasproperty(trans,:t)
    base=affine ? trans.t : trans
    nameof(typeof(base))===:Agnesi2Transform || error("unsupported source transform")
    poly=radial.J; env=radial.envelope
    nameof(typeof(env)) in (:OneEnvelope,:PolyEnvelope) || error("unsupported source envelope")
    (;a=base.a,rin=base.rin,r0=base.r0,p=Int(base.p),q=Int(base.q),
      x1=affine ? trans.x1 : 0.,x2=affine ? trans.x2 : 1.,y1=affine ? trans.y1 : 0.,y2=affine ? trans.y2 : 1.,
      A=copy(poly.A),B=copy(poly.B),C=copy(poly.C),pl=poly.pl,pr=poly.pr,tl=poly.tl,tr=poly.tr,
      envelope_power=hasproperty(env,:p) ? env.p : 0,envelope_r0=hasproperty(env,:r0) ? env.r0 : 1.,offset,cutoff)
end
function channel(radial,z,z0,offset,nr,cutoff)
    c=config(radial,z,z0,offset,cutoff)
    metadata=Dict("spec"=>[[n,0] for n=1:nr],"cutoff"=>cutoff,"index_weights"=>Int[],"weights_known"=>false,
      "simd"=>false,"ad"=>false,"origin"=>"error","envelope"=>"included","channel_dependence"=>"central_neighbor")
    refs=Dict[]; ev=0.; ed=0.
    tmp=J.alloc_temp(radial); tmpd=J.alloc_temp_d(radial)
    for r in range(.45,cutoff-.0001;length=137)
        v=zeros(nr); d=similar(v); LegacyPolynomial()(v,d,r,NamedTuple(),c)
        original=zeros(length(radial)); derivative=zeros(length(radial))
        J.evaluate_d!(original,derivative,tmpd,radial,r,z,z0)
        J.evaluate!(original,tmp,radial,r,z,z0)
        ev=max(ev,maximum(abs,v[offset+1:offset+length(original)]-original))
        ed=max(ed,maximum(abs,d[offset+1:offset+length(original)]-derivative))
    end
    ev<1e-9 && ed<1e-7 || error("radial conversion mismatch $ev / $ed")
    for fraction in (.23,.59,.83)
        r=fraction*cutoff; v=zeros(nr); d=similar(v); LegacyPolynomial()(v,d,r,NamedTuple(),c)
        push!(refs,Dict("r"=>r,"values"=>v,"derivatives"=>d))
    end
    Dict("kind"=>"callable","metadata"=>metadata,"config"=>encode(c),"route"=>"registered","interface_version"=>1,
         "id"=>"PotentialFrameWork.legacy_polynomial.v1","reference"=>refs),ev,ed
end
function main(source,destination)
    raw=JSON.parsefile(source); ip=ACEpotentials.load_potential(source)
    length(ip.components)==3 || error("expected pair, PI, reference")
    pair,pi,onebody=ip.components
    pb=pi.pibasis; b1=pb.basis1p
    nameof(typeof(b1))===:BasicPSH1pBasis || error("unsupported source one-particle basis")
    z=collect(pb.zlist.list); nz=length(z); nrad=length(b1.J); npair=length(pair.basis.J[1])
    lmax=maximum(b.l for b in b1.spec); cutoff=maximum(Float64(J.cutoff(v)) for v in (pair,pi))
    densities=Tuple{Int,Int}[]; density_ids=Dict{Tuple{Int,Int},Int}()
    products=Vector{Int}[]; product_ids=Dict{Tuple,Int}()
    rows=Int[]; cols=Int[]; vals=Float64[]; coeffs=Vector{Float64}[]
    # Each source PI coefficient stays explicit and in original central-species
    # order. Real(C-product) is expanded exactly into real harmonics; no fitting.
    nrows=sum(length,pi.coeffs); WB=zeros(nrows,nz); row=0
    for iz=1:nz
        specs=Old.get_basis_spec(pb,iz)
        length(specs)==length(pi.coeffs[iz]) || error("source PI coordinate mismatch")
        for (j,spec) in enumerate(specs)
            row+=1; WB[row,iz]=pi.coeffs[iz][j]
            terms=Dict{Tuple,ComplexF64}(() => 1.0 + 0im)
            for b in spec.oneps
                jz=Old.z2i(pb,b.z); radial=b.n+(jz-1)*nrad; l=b.l; m=b.m
                choices=m==0 ? [(l*l+l+1,1.0 + 0im)] :
                  [(l*l+l+abs(m)+1,(m>0 ? (-1.0)^m : 1.0)/sqrt(2)+0im),
                   (l*l+l-abs(m)+1,(m>0 ? (-1.0)^m : -1.)*im/sqrt(2))]
                next=Dict{Tuple,ComplexF64}()
                for (old,w) in terms,(angular,factor) in choices
                    key=(radial,angular)
                    a=get!(density_ids,key) do; push!(densities,key); length(densities); end
                    factors=Tuple(sort!([old...,a])); next[factors]=get(next,factors,0im)+w*factor
                end
                terms=next
            end
            for (indices,w) in terms
                iszero(real(w)) && continue
                isempty(indices) && error("constant source PI term requires explicit reference mapping")
                k=get!(product_ids,indices) do; push!(products,collect(indices)); length(products); end
                push!(rows,row); push!(cols,k); push!(vals,real(w))
            end
        end
    end
    radials=Dict[]; pairs=Dict[]; ev=0.; ed=0.; Wpair=zeros(npair*nz,nz); pair_map=zeros(Int,npair*nz,nz)
    for jz=1:nz,iz=1:nz
        ch,v,d=channel(b1.J,z[jz],z[iz],(jz-1)*nrad,nrad*nz,cutoff); push!(radials,ch); ev=max(ev,v); ed=max(ed,d)
        ch,v,d=channel(pair.basis.J[jz,iz],z[jz],z[iz],(jz-1)*npair,npair*nz,cutoff); push!(pairs,ch); ev=max(ev,v); ed=max(ed,d)
        first=pair.basis.bidx0[jz,iz]
        Wpair[(jz-1)*npair+1:jz*npair,iz]=pair.coeffs[first+1:first+npair]/2
        pair_map[(jz-1)*npair+1:jz*npair,iz]=collect(first+1:first+npair)
    end
    family(ch,n)=Dict("kind"=>"interface_channels_v1","spec"=>[Dict("n"=>i,"l"=>0) for i=1:n],"shape"=>[nz,nz],"channels"=>ch)
    resolved=Dict("schema"=>"PotentialFrameWork.resolved_ACE","schema_version"=>1,"source_sha256"=>bytes2hex(sha256(read(source))),
      "species_atomic_numbers"=>Int.(z),"reference_eV"=>[onebody.E0[J.chemical_symbol(zz)] for zz in z],
      "angular"=>Dict("kind"=>"real_spherical","normalization"=>"L2","ordering"=>"l*l+l+m+1","lmax"=>lmax),
      "radial"=>family(radials,nrad*nz),"pair_radial"=>family(pairs,npair*nz),"density_indices"=>collect.(densities),"products"=>products,
      "coupling"=>Dict("shape"=>[nrows,length(products)],"rows"=>rows,"columns"=>cols,"values"=>vals),
      "parameters"=>Dict("WB_shape"=>collect(size(WB)),"Wpair_shape"=>collect(size(Wpair)),"flat"=>vcat(vec(WB),vec(Wpair)),
          "ordering"=>"vec(WB), vec(Wpair); Julia column-major"),
      "conversion"=>Dict("kind"=>"exact_legacy_complex_PI_to_real","radial_value_max_error"=>ev,"radial_derivative_max_error"=>ed,
         "producer_julia"=>string(VERSION),"ACEpotentials"=>string(pkgversion(ACEpotentials)),"original_PI_per_species"=>collect(length.(pi.coeffs)),"original_PI_species_offsets"=>vcat(0,cumsum(collect(length.(pi.coeffs)))),
         "Wpair_source_indices_column_major"=>vec(pair_map),"Wpair_source_multiplier"=>0.5,
         "limitations"=>"Source JSON contains precontracted PI coordinates, not the original invariant fitting basis. Individual rows are not rotation invariants; invariant-row selection/fitting on this representation is unsupported."))
    open(destination,"w") do io; TOML.print(io,resolved;sorted=true); end
    println(destination," densities=",length(densities)," products=",length(products)," PI coordinates=",nrows," radial errors=",(ev,ed))
end
length(ARGS)==2 || error("usage: export_legacy.jl SOURCE.json OUTPUT.toml")
main(ARGS...)
