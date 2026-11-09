# Numerical evidence on the actual AMD device, independent energy stencils.
using PotentialFrameWork, AMDGPU, LinearAlgebra, TOML
const A=PotentialFrameWork.ACE
const AB=A.AtomsBase
const U=A.Unitful
AMDGPU.functional() || error("AMD validation requested but unavailable")
AMDGPU.device!(only(filter(d->occursin("6800 XT",string(d)),AMDGPU.devices())))
model,ps,_=A.load_resolved(joinpath(@__DIR__,"../../data/ace_ref/resolved_model_v1.toml"))
bundle=A.adapt_model(model,ps;array_type=AMDGPU.ROCArray)
function system(X,H,pb)
    cells=ntuple(j->A.SVector{3}(H[:,j])*U.u"Å",3)
    AB.atomic_system([[:Ti,:Al,:Ti][i]=>A.SVector{3}(X[:,i])*U.u"Å" for i=1:3],cells,pb)
end
X=[.4 2.6 1.1;0. .2 2.8;0. .3 .1]; H=Matrix{Float64}(I,3,3)*18
sys=system(X,H,(false,false,false))
small_calculator=A.ACECalculator(model,ps;array_type=AMDGPU.ROCArray,precision=Float32)
@assert typeof(A.AC.zero_energy(sys,small_calculator).val)===Float32
executor=A.GPUExecutor(sys,bundle.model,bundle.parameters)
y=A.gpu_result(A.execute_gpu!(executor)); h=1e-5
energy(s)=only(A.gpu_result(A.execute_gpu!(A.prepare_executor!(executor,s);properties=(A.AC.Energy(),))).energy)
force_error=0.; strain_error=0.
for a=1:3,i=1:3
    plus=copy(X); minus=copy(X); plus[a,i]+=h; minus[a,i]-=h
    fd=-(energy(system(plus,H,(false,false,false)))-energy(system(minus,H,(false,false,false))))/(2h)
    global force_error=max(force_error,abs(fd-y.forces[a,i]))
end
for a=1:3,b=1:3
    delta=zeros(3,3); delta[b,a]=h
    fd=-(energy(system((I+delta)*X,H,(false,false,false)))-energy(system((I-delta)*X,H,(false,false,false))))/(2h)
    global strain_error=max(strain_error,abs(fd-y.virial[a,b]))
end
periodicX=copy(X); periodicX[1,1]=17.8
periodic=system(periodicX,H,(true,true,true)); A.prepare_executor!(executor,periodic)
p=A.gpu_result(A.execute_gpu!(executor)); periodic_strain_error=0.
for a=1:3
    delta=zeros(3,3); delta[a,a]=h
    fd=-(energy(system((I+delta)*periodicX,(I+delta)*H,(true,true,true)))-energy(system((I-delta)*periodicX,(I-delta)*H,(true,true,true))))/(2h)
    global periodic_strain_error=max(periodic_strain_error,abs(fd-p.virial[a,a]))
end
projected=A.project_species(model,[:Al],0:3)
projected_bundle=A.adapt_model(projected,ps;array_type=AMDGPU.ROCArray)
projected_y=A.gpu_result(A.execute_gpu!(A.GPUExecutor(sys,projected_bundle.model,projected_bundle.parameters)))
projection_force_error=maximum(abs,projected_y.forces-y.forces)
report=Dict("julia"=>string(VERSION),"device"=>string(AMDGPU.device()),"step"=>h,"absolute_tolerance"=>2e-7,
    "force_fd_max_eV_per_A"=>force_error,"open_strain_fd_max_eV"=>strain_error,"periodic_diagonal_strain_fd_max_eV"=>periodic_strain_error,
    "net_force_norm_eV_per_A"=>norm(vec(sum(y.forces;dims=2))),"full_order_projection_force_max_eV_per_A"=>projection_force_error)
path=isempty(ARGS) ? "derivative-checks.toml" : ARGS[1]
open(path,"w") do io; TOML.print(io,report;sorted=true); end
println(report)

@assert force_error<=2e-7 && strain_error<=2e-7 && periodic_strain_error<=2e-7
@assert projection_force_error<=1e-9
