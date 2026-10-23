module PotentialFrameWorkSpheriCartExt
using PotentialFrameWork, SpheriCart
import PotentialFrameWork.ACE: angular_basis, angular_values, angular_ed, angular_ed!
import PotentialFrameWork.ACE: angular_basis_generated
struct ConstantAngular{T}
    value::T
end
angular_basis(lmax,::Type{T}) where T = lmax==0 ? ConstantAngular(T(inv(sqrt(4pi)))) : SpheriCart.SphericalHarmonics(lmax;normalisation=:L2,static=false,T)
struct GeneratedAngular{B}
    basis::B
end
angular_basis_generated(lmax,::Type{T}) where T = lmax==0 ? ConstantAngular(T(inv(sqrt(4pi)))) : GeneratedAngular(SpheriCart.SphericalHarmonics(lmax;normalisation=:L2,static=true,T))
function angular_ed!(Y,dY,b::ConstantAngular,Rs)
    fill!(Y,b.value); fill!(dY,zero(eltype(dY)))
    Y,dY
end
function angular_ed(b::ConstantAngular,Rs)
    T=eltype(eltype(Rs))
    fill(b.value+zero(T),length(Rs),1),fill(zero(eltype(Rs)),length(Rs),1)
end
function angular_ed!(Y,dY,b::GeneratedAngular,Rs)
    for i in eachindex(Rs)
        values,gradients=SpheriCart.compute_with_gradients(b.basis,Rs[i])
        for n in eachindex(values); Y[i,n]=values[n]; dY[i,n]=gradients[n]; end
    end
    Y,dY
end
angular_values(basis::SpheriCart.SphericalHarmonics,rs)=SpheriCart.compute(basis,rs)
angular_ed(basis::SpheriCart.SphericalHarmonics,rs)=SpheriCart.compute_with_gradients(basis,rs)
angular_ed!(Y,dY,basis::SpheriCart.SphericalHarmonics,rs)=SpheriCart.compute_with_gradients!(Y,dY,basis,rs)
end
