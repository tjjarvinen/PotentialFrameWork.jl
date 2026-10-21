module Execution
using Unitful
using ..NeighbourList
using ..NeighbourList: lane
using ..PairPotentials
using ..PairPotentials: select
export evaluate, evaluate!, execute!, outputs

@inline function displacement(d,L,periodic)
    periodic ? d-round(d/L)*L : d
end

@inline function geometry(w::Workspace{W,T}, i, firstj, lastj, image=nothing) where {W,T}
    js=ntuple(k -> min(firstj+k-1,lastj),Val(W)) # every index valid, including padding
    self = image !== nothing && image.self
    structural=pack(Val(W),k -> firstj+k-1<=lastj && (js[k]>i || (self && js[k]==i)))
    shift = image === nothing ? (zero(T),zero(T),zero(T)) : image.shift
    dx=pack(Val(W), k -> (w.minimum_image ? displacement(w.x[js[k]]-w.x[i],w.lengths[1],w.periodic[1]) : w.x[js[k]]-w.x[i]+shift[1]))
    dy=pack(Val(W), k -> (w.minimum_image ? displacement(w.y[js[k]]-w.y[i],w.lengths[2],w.periodic[2]) : w.y[js[k]]-w.y[i]+shift[2]))
    dz=pack(Val(W), k -> (w.minimum_image ? displacement(w.z[js[k]]-w.z[i],w.lengths[3],w.periodic[3]) : w.z[js[k]]-w.z[i]+shift[3]))
    active=structural & (dx*dx+dy*dy+dz*dz < w.cutoff*w.cutoff)
    # Safe geometry BEFORE callbacks (including divisions/logs), not NaN*0 after.
    dx=select(active,dx,one(dx)); dy=select(active,dy,zero(dy)); dz=select(active,dz,zero(dz))
    block=PairBlock(w.species[i],ntuple(k -> w.species[js[k]],Val(W)),dx,dy,dz,structural,active)
    js,block
end

@inline function update_virial(vir, d, f)
    ntuple(k -> vir[k]-sum(d[mod1(k,3)]*f[cld(k,3)]),Val(9))
end

@inline function cell_task!(w::Workspace{W,T}, p, req, task) where {W,T}
    e=zero(T)
    vir=ntuple(_ -> zero(T),9)
    for i in w.ranges[task]
        central=ntuple(_ -> zero(T),3)
        for image in w.images[task]
            r=w.ranges[image.cell]
            # Same-cell half-list skips the already owned lower triangle entirely.
            lo=max(first(r),i+(image.self ? 0 : 1))
            for j=lo:W:last(r)
                js,b=geometry(w,i,j,last(r),image)
                any(b.active) || continue
                out=p(b,req)
                en=select(b.active,out.energy,zero(out.energy))
                e+=sum(en)
                if req isa EnergyForcesVirial
                    fx=select(b.active,out.fx,zero(out.fx))
                    fy=select(b.active,out.fy,zero(out.fy))
                    fz=select(b.active,out.fz,zero(out.fz))
                    central=(central[1]+sum(fx),central[2]+sum(fy),central[3]+sum(fz))
                    # Scalar endpoint updates deliberately handle duplicate padding destinations.
                    for k=1:W
                        lane(b.active,k) || continue
                        dest=js[k]
                        w.force[1,dest]-=lane(fx,k)
                        w.force[2,dest]-=lane(fy,k)
                        w.force[3,dest]-=lane(fz,k)
                    end
                    f=(fx,fy,fz); d=(b.dx,b.dy,b.dz)
                    vir=update_virial(vir,d,f)
                end
            end
        end
        if req isa EnergyForcesVirial
            for a=1:3; w.force[a,i]+=central[a]; end
        end
    end
    w.energies[task]=e
    for a=1:9; w.virials[a,task]=vir[a]; end
    nothing
end

"""Execute already prepared numeric storage. Internal benchmarking entrypoint;
positions/species/cell must not change between prepare! and this call.
Each phase is a barrier; each cell task is serial and owns its scalar partials.
"""
function execute!(w, bound, req; parallel=false)
    w.schedule_version==w.version || error("schedule topology version mismatch")
    fill!(w.force,0); fill!(w.energies,0); fill!(w.virials,0)
    if parallel
        for phase in w.phases
            Threads.@threads :dynamic for k in eachindex(phase)
                cell_task!(w,bound,req,phase[k])
            end
        end
    else
        for task in eachindex(w.ranges); cell_task!(w,bound,req,task); end
    end
    w
end

function outputs(w::Workspace{W,T},req) where {W,T}
    t=time_ns()
    energy=sum(w.energies)*u"eV"
    if req isa EnergyOnly
        w.timings[:partial_merge]=(time_ns()-t)*1e-9
        w.timings[:output_permutation]=0.0
        return (;energy)
    end
    v=reshape(vec(sum(w.virials;dims=2)),3,3)*u"eV"
    w.timings[:partial_merge]=(time_ns()-t)*1e-9
    t=time_ns()
    f=Matrix{typeof(zero(T)*u"eV/Å")}(undef,3,length(w.p))
    for k in eachindex(w.p), a=1:3
        f[a,w.p[k]]=w.force[a,k]*u"eV/Å"
    end
    w.timings[:output_permutation]=(time_ns()-t)*1e-9
    (;energy,forces=f,virial=v)
end

"""Evaluate an AtomsBase system, rebuilding all geometry/topology and returning
unitful original-order outputs. `evaluate!` reuses workspace storage, never the
caller's positions. Do not share one mutable workspace between concurrent calls.
"""
function evaluate!(w::Workspace,sys,p; request=EnergyForcesVirial(), parallel=false)
    request isa Union{EnergyOnly,EnergyForcesVirial} || throw(ArgumentError("unsupported request"))
    prepare!(w,sys,cutoff(p))
    t=time_ns()
    bound=bind_species(p,w.species_map)
    w.timings[:parameters]=(time_ns()-t)*1e-9
    t=time_ns()
    execute!(w,bound,request;parallel)
    w.timings[:evaluation]=(time_ns()-t)*1e-9
    result=outputs(w,request)
    result
end
function evaluate(sys,p; width=4, precision=Float64, bruteforce_threshold=32, kwargs...)
    evaluate!(Workspace(;width,precision,bruteforce_threshold),sys,p;kwargs...)
end
end
