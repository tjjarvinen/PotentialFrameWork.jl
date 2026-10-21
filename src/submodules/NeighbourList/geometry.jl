# Scalar/SIMD masking is shared by neighbour packing and pair callbacks.
@inline select(mask::Bool,a,b) = ifelse(mask,a,b)
@inline select(mask::Vec,a,b) = vifelse(mask,a,b)

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

