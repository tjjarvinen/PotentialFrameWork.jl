module PotentialFrameWorkForwardDiffExt
using PotentialFrameWork, ForwardDiff
import PotentialFrameWork: energy_parameter_gradient
import PotentialFrameWork.PairPotentials: pair_vjp
using PotentialFrameWork.PairPotentials: edge_result, parameter_seed
using StaticArrays: SVector
using LinearAlgebra: dot
using PotentialFrameWork.ACE: ResolvedModel, site_data, site_reverse, site_parameter_gradient, edge_displacement
using PotentialFrameWork.PairPotentials: parameter_flatten, parameter_rebuild
import PotentialFrameWork.ACE: ace_vjp, custom_parameter_gradient, radial_ad!, radial_ad_probe

# Preserve nested field names and scalar/array shapes at the AD boundary.
unpack(template::NamedTuple,flat)=parameter_rebuild(template,flat)
function energy_parameter_gradient(sys,p::AbstractPairPotential,parameters::NamedTuple;kwargs...)
    PotentialFrameWork.PairPotentials.parameter_seed(parameters)
    flat=parameter_flatten(parameters)
    isempty(flat) && return unpack(parameters,flat)
    gradient=ForwardDiff.gradient(flat) do theta
        parameter_energy(sys,p,unpack(parameters,theta);kwargs...)
    end
    unpack(parameters,gradient)
end

# Edge-local AD is a correctness fallback: storage is O(N + parameters), with
# parameter work proportional to edges × parameter count. No system Jacobian or
# Hessian is built. Custom callback force components must be differentiable.
function pair_vjp(X,p,ps,env,H,de,df,dw)
    flat=parameter_flatten(ps)
    dp=zero(flat)
    dx=zeros(promote_type(eltype(X),eltype(flat)),size(X))
    bound=bind_species(with_parameters(p,ps),unique(env.labels))
    for k in eachindex(env.centres)
        i=env.centres[k]; j=env.neighbours[k]
        d=SVector{3}(ntuple(a -> X[a,j]-X[a,i]+sum(H[a,b]*env.shifts[k][b] for b=1:3),3))
        sum(abs2,d)<env.radius^2 || continue
        b=SVector{3}(ntuple(a -> df[a,i]-df[a,j],3))
        function loss(out,r)
            g=(out.fx,out.fy,out.fz)
            de*out.energy+sum(b[a]*g[a] for a=1:3)-sum(dw[a,c]*r[a]*g[c] for a=1:3,c=1:3)
        end
        gd=ForwardDiff.gradient(r -> loss(edge_result(bound,env,k,r,EnergyForcesVirial()),r),d)
        gp=ForwardDiff.gradient(flat) do theta
            params=unpack(ps,theta)
            updated=bind_species(with_parameters(p,params),unique(env.labels))
            promoted=d .+ parameter_seed(params)
            loss(edge_result(updated,env,k,promoted,EnergyForcesVirial()),promoted)
        end
        dp .+= gp
        for a=1:3
            dx[a,i]-=gd[a]; dx[a,j]+=gd[a]
        end
    end
    (;positions=dx,parameters=unpack(ps,dp))
end

struct ACEParameterDirection end
# The saved displacement direction is the outer operation's tangent. Local
# parameter/distance differentiation nests outside it and must retain that
# direction in the resulting mixed product. Only our own tag is ordered here.
ForwardDiff.:(≺)(::Type{ACEParameterDirection},::Type{<:ForwardDiff.Tag})=true
ForwardDiff.:(≺)(::Type{<:ForwardDiff.Tag},::Type{ACEParameterDirection})=false

function ace_vjp_impl(X,model,ps,env,H,de,df,dw;centers=nothing)
    for family in (model.radial,model.pair_radial),ch in family; PotentialFrameWork.ACE.require_radial_ad(ch); end
    flat=parameter_flatten(ps); dp=zero(flat); dx=zeros(eltype(X),size(X))
    basis=PotentialFrameWork.ACE.angular_basis(model.lmax,eltype(model.reference))
    for i in PotentialFrameWork.ACE.center_indices(centers,length(env.species))
        iz=env.species[i]
        edges=[k for k=env.offsets[1,i]:env.offsets[end,i]-1 if env.active[k]]
        isempty(edges) && continue
        Rs=[edge_displacement(X,H,env,i,k) for k in edges]
        zs=[env.species[env.neighbours[k]] for k in edges]
        bs=[SVector{3}(df[:,i]-df[:,env.neighbours[k]]) for k in edges]
        # For v_j = B_i-B_j-C'R_j, the displacement adjoint is
        # de*g + Hessian(E)*v - C*g. Hessian symmetry allows one directional
        # differentiation of the analytic reverse, independent of neighbour or
        # parameter count. The final term differentiates virial geometry.
        dirs=[bs[j]-transpose(dw)*Rs[j] for j in eachindex(Rs)]
        dualrs=[SVector{3}(ntuple(a -> ForwardDiff.Dual{ACEParameterDirection}(Rs[j][a],dirs[j][a]),3)) for j in eachindex(Rs)]
        data=site_data(model,ps,iz,zs,dualrs;basis)
        reverse=site_reverse(model,ps,iz,data)
        for (j,k) in enumerate(edges)
            gvalue=ForwardDiff.value.(reverse.gradients[j])
            action=map(g -> ForwardDiff.partials(g)[1],reverse.gradients[j])
            gd=de*gvalue+action-dw*gvalue
            for a=1:3
                dx[a,i]-=gd[a]; dx[a,env.neighbours[k]]+=gd[a]
            end
        end
        # The same directional pass supplies the mixed parameter product.
        g=parameter_flatten(site_parameter_gradient(model,ps,iz,zs,data;reverse))
        for k in eachindex(dp)
            dp[k]+=de*ForwardDiff.value(g[k])+ForwardDiff.partials(g[k])[1]
        end
    end
    (;positions=dx,parameters=parameter_rebuild(ps,dp))
end
function ace_vjp(X,model,ps,env,H,de,df,dw;kwargs...)
    PotentialFrameWork.ACE.has_source(model) ? Base.invokelatest(ace_vjp_impl,X,model,ps,env,H,de,df,dw;kwargs...) : ace_vjp_impl(X,model,ps,env,H,de,df,dw;kwargs...)
end

# Scalar values-only adapter. This correctness fallback allocates temporary
# channel vectors per derivative sample; it is not the allocation-free SIMD API.
function radial_ad!(v,d,f,r::Real,ps,c)
    f(v,r,ps,c)
    for n in eachindex(v)
        d[n]=ForwardDiff.derivative(r) do x
            tmp=similar(v,typeof(x+PotentialFrameWork.ACE.radial_parameter_seed(ps)))
            f(tmp,x,ps,c); tmp[n]
        end
    end
    nothing
end
function custom_parameter_gradient(model,ps,iz,zs,data,g::NamedTuple)
    result=g
    for key in (:custom_radial,:custom_pair)
        hasproperty(ps,key) || continue
        template=getproperty(ps,key); flat=parameter_flatten(template)
        isempty(flat) && continue
        gradient=ForwardDiff.gradient(flat) do values
            params=merge(ps,NamedTuple{(key,)}((parameter_rebuild(template,values),)))
            site_data(model,params,iz,zs,data.Rs).energy
        end
        result=merge(result,NamedTuple{(key,)}((parameter_rebuild(template,gradient),)))
    end
    result
end

function radial_ad_probe(ch::PotentialFrameWork.ACE.CallableRadial,ps,r)
    A=PotentialFrameWork.ACE
    try
        _,d=A.radial_ed(ch,r,ps)
        first_derivative=ForwardDiff.jacobian(x -> A.radial_values(ch,x[1],ps),[r])[:,1]
        second=ForwardDiff.jacobian(x -> last(A.radial_ed(ch,x[1],ps)),[r])
        flat=parameter_flatten(ps)
        mixed=isempty(flat) ? Float64[] : ForwardDiff.gradient(p -> sum(last(A.radial_ed(ch,r,parameter_rebuild(ps,p)))),flat)
        passed=isapprox(d,first_derivative;atol=1e-9,rtol=1e-7) && all(isfinite,second) && all(isfinite,mixed)
        return (;checked=true,passed,reason=passed ? "" : "local radial AD derivative consistency failed")
    catch e
        return (checked=true,passed=false,reason="local radial AD: "*sprint(showerror,e))
    end
end

end
