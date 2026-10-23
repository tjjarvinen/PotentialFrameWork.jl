linear_parameters(ps)=vcat(vec(ps.WB),vec(ps.Wpair))
function install_readout(ps,coefficients;fit_pair=true)
    n=length(ps.WB)
    if !fit_pair
        length(coefficients)==n || throw(DimensionMismatch("ACE readout length"))
        return merge(deepcopy(ps),(;WB=reshape(copy(coefficients),size(ps.WB))))
    end
    length(coefficients)==n+length(ps.Wpair) || throw(DimensionMismatch("readout length"))
    merge(ps,(;WB=reshape(copy(coefficients[1:n]),size(ps.WB)),
               Wpair=reshape(copy(coefficients[n+1:end]),size(ps.Wpair))))
end

"""Dense E/F/W linear design on fixed internal radial parameters. Rows are
energy, `vec(forces)`, `vec(virial)` in original atom order, with Julia column
ordering. `baseline` must be subtracted from labels before fitting. Imported
pair columns and ACE columns are jointly linear in the original coordinates.
"""
function design_matrix(X,model,ps,env,H;fit_pair=true)
    nb=size(model.coupling,1); nq=length(model.pair_spec); nz=length(model.atomic_numbers)
    n=length(env.species); cols=(nb+nq)*nz
    A=zeros(promote_type(eltype(X),eltype(ps.WB)),1+3n+9,cols)
    baseline=zeros(eltype(A),size(A,1)); baseline[1]=sum(model.reference[env.species])
    basis=angular_basis(model.lmax,eltype(model.reference))
    for i=1:n
        iz=env.species[i]; data=site_diagnostics(X,model,ps,env,H,i;basis)
        wb=(iz-1)*nb .+ (1:nb); wp=nb*nz+(iz-1)*nq .+ (1:nq)
        A[1,wb] .+= data.B
        A[1,wp] .+= vec(sum(data.Q;dims=1))
        # Differentiate products once in density coordinates, then apply the
        # fixed sparse invariant map. Never evaluate the whole model per column.
        D=zeros(eltype(A),length(model.products),length(model.density))
        for (k,indices) in enumerate(model.products), t in eachindex(indices)
            v=one(eltype(A))
            for u in eachindex(indices); u==t || (v*=data.A[indices[u]]); end
            D[k,indices[t]]+=v
        end
        dBdA=model.coupling*D
        for (j,k) in enumerate(data.edges)
            d=data.Rs[j]; direction=d/data.rs[j]
            dAdd=zeros(eltype(A),length(model.density),3)
            for (a,(r,m)) in enumerate(model.density), axis=1:3
                dAdd[a,axis]=data.dR[j,r]*data.Y[j,m]*direction[axis]+data.R[j,r]*data.dY[j,m][axis]
            end
            dB=dBdA*dAdd
            dQ=data.dQ[j,:]*transpose(direction)
            dest=env.neighbours[k]
            for (indices,derivatives) in ((wb,dB),(wp,dQ)), a=1:3
                A[1+3(i-1)+a,indices] .+= derivatives[:,a]
                A[1+3(dest-1)+a,indices] .-= derivatives[:,a]
                for b=1:3
                    row=1+3n+a+3(b-1)
                    A[row,indices] .-= d[a]*derivatives[:,b]
                end
            end
        end
    end
    if !fit_pair
        baseline .+= A[:,nb*nz+1:end]*vec(ps.Wpair)
        return (;matrix=A[:,1:nb*nz],baseline)
    end
    (;matrix=A,baseline)
end

"""Diagonal algebraic prior in original invariant/readout coordinates.
Each invariant uses the largest sum(n+l) among its structurally coupled products;
pair columns use (1+n+l)^p. No precontracted coordinate prior is substituted.
"""
function smoothness_prior(model;p=4,fit_pair=true)
    degrees=[sum(sum(model.radial_spec[model.density[a][1]]) for a in product) for product in model.products]
    d=ones(Float64,size(model.coupling,1))
    for j=1:size(model.coupling,2), k in nzrange(model.coupling,j)
        i=rowvals(model.coupling)[k]; d[i]=max(d[i],1+degrees[j])
    end
    nz=length(model.atomic_numbers)
    Diagonal(fit_pair ? vcat(repeat(d.^p,nz),repeat([(1+sum(s))^p for s in model.pair_spec],nz)) : repeat(d.^p,nz))
end

function fit_linear(args...;kwargs...)
    throw(ArgumentError("load ACEfit for linear regression"))
end

"""Select available numeric eV/Å labels from a design result; a missing property
or `nothing` omits those rows. Property weights multiply residuals directly.
Returns baseline-subtracted `labels` and matching row weights for `fit_linear`.
Force labels are 3×N and virial labels are 3×3, in original atom order.
"""
function training_rows(design,labels;weights=(energy=1.,forces=1.,virial=1.))
    n3=size(design.matrix,1)-10
    indices=Int[]; values=eltype(design.matrix)[]; residual_weights=Float64[]
    for (key,range) in ((:energy,1:1),(:forces,2:1+n3),(:virial,2+n3:10+n3))
        hasproperty(labels,key) || continue
        value=getproperty(labels,key); isnothing(value) && continue
        v=key===:energy ? [value] : vec(value)
        length(v)==length(range) || throw(DimensionMismatch("$key label dimensions"))
        weight=getproperty(weights,key)
        isfinite(weight) && weight>=0 || throw(ArgumentError("invalid $key weight"))
        append!(indices,range); append!(values,v); append!(residual_weights,fill(weight,length(range)))
    end
    (;matrix=design.matrix[indices,:],labels=values-design.baseline[indices],row_weights=residual_weights)
end
