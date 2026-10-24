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
function design_matrix_impl(X,model,ps,env,H;fit_pair=true,centers=nothing)
    nb=size(model.coupling,1); nq=length(model.pair_spec); nz=length(model.atomic_numbers)
    n=length(env.species); cols=(nb+nq)*nz
    A=zeros(promote_type(eltype(X),eltype(ps.WB)),1+3n+9,cols)
    baseline=zeros(eltype(A),size(A,1)); selected_centers=center_indices(centers,n)
    baseline[1]=sum((reference_value(model,env.species[i]) for i in selected_centers);init=zero(eltype(A)))
    basis=angular_basis(model.lmax,eltype(model.reference))
    for i in selected_centers
        iz=env.species[i]; data=site_diagnostics_impl(X,model,ps,env,H,i;basis)
        wb=(iz-1)*nb .+ (1:nb); wp=nb*nz+(iz-1)*nq .+ (1:nq)
        A[1,wb] .+= data.B
        A[1,wp] .+= vec(sum(data.Q;dims=1))
        # Differentiate products once in density coordinates, then apply the
        # fixed sparse invariant map. Never evaluate the whole model per column.
        dBdU=zeros(eltype(A),nb,length(model.density)); dBdV=similar(dBdU)
        for row=1:nb
            weights=Vector(model.coupling[row,:]).*model.policy.mask[row]
            proj=density_reverse(model,weights,data.U,data.V)
            dBdU[row,:]=proj.dU; dBdV[row,:]=proj.dV
        end
        for (j,k) in enumerate(data.edges)
            d=data.Rs[j]; direction=d/data.rs[j]
            dAdd=zeros(eltype(A),length(model.density),3)
            for (a,(r,m)) in enumerate(model.density), axis=1:3
                dAdd[a,axis]=data.dR[j,r]*data.Y[j,m]*direction[axis]+data.R[j,r]*data.dY[j,m][axis]
            end
            dB=(data.zs[j] in model.policy.species ? dBdV : dBdU)*dAdd
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
function smoothness_prior(model;p=4,fit_pair=true,unknown=:error)
    radial_weights=family_index_weights(model)
    pair_weights=fit_pair ? family_index_weights(model;pair=true) : Int[]
    if isnothing(radial_weights) || isnothing(pair_weights)
        unknown===:neutral || throw(ArgumentError("undefined custom index weights: supply a prior or unknown=:neutral"))
        n=(size(model.coupling,1)+(fit_pair ? length(model.pair_spec) : 0))*length(model.atomic_numbers)
        return Diagonal(ones(n))
    end
    degrees=[sum(radial_weights[model.density[a][1]] for a in product) for product in model.products]
    d=ones(Float64,size(model.coupling,1))
    for j=1:size(model.coupling,2), k in nzrange(model.coupling,j)
        i=rowvals(model.coupling)[k]; d[i]=max(d[i],1+degrees[j])
    end
    nz=length(model.atomic_numbers)
    Diagonal(fit_pair ? vcat(repeat(d.^p,nz),repeat((1 .+ pair_weights).^p,nz)) : repeat(d.^p,nz))
end

function fit_linear(args...;kwargs...)
    throw(ArgumentError("load ACEfit for linear regression"))
end

training_label_value(value,unit)=value
training_label_value(value::Unitful.AbstractQuantity,unit)=ustrip(unit,value)

"""Select available labels from a design result; a missing property
or `nothing` omits those rows. Property weights multiply residuals directly.
Returns baseline-subtracted `labels` and matching row weights for `fit_linear`.
Unitful labels are converted to eV (energy/virial) and eV/Å (forces) before
baseline subtraction. Plain numbers use those same units. Incompatible dimensions
are rejected. Weights apply to residuals in these converted units.
Force labels are 3×N or a vector of N three-component vectors, in original atom
order; virial labels are 3×3 configurational virials (energy, without volume
normalization). Geometry is converted separately by `ace_inputs`.
"""
function training_rows(design,labels;weights=(energy=1.,forces=1.,virial=1.))
    n3=size(design.matrix,1)-10
    indices=Int[]; values=eltype(design.matrix)[]; residual_weights=Float64[]
    for (key,range) in ((:energy,1:1),(:forces,2:1+n3),(:virial,2+n3:10+n3))
        hasproperty(labels,key) || continue
        value=getproperty(labels,key); isnothing(value) && continue
        v=if key===:energy
            (value,)
        elseif key===:forces && value isa AbstractVector && all(x -> x isa AbstractVector,value)
            all(x -> length(x)==3,value) || throw(DimensionMismatch("forces need three components per atom"))
            collect(Iterators.flatten(value))
        else
            vec(value)
        end
        length(v)==length(range) || throw(DimensionMismatch("$key label dimensions"))
        weight=getproperty(weights,key)
        isfinite(weight) && weight>=0 || throw(ArgumentError("invalid $key weight"))
        unit=key===:forces ? u"eV/Å" : u"eV"
        append!(indices,range)
        append!(values,(training_label_value(x,unit) for x in v))
        append!(residual_weights,fill(weight,length(range)))
    end
    (;matrix=design.matrix[indices,:],labels=values-design.baseline[indices],row_weights=residual_weights)
end

function design_matrix(X,model,ps,env,H;kwargs...)
    has_source(model) ? Base.invokelatest(design_matrix_impl,X,model,ps,env,H;kwargs...) : design_matrix_impl(X,model,ps,env,H;kwargs...)
end
