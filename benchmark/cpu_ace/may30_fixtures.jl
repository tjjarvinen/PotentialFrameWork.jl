# Pure deterministic geometry, shared by the live comparison and core regressions.
# Columns are lattice generators; all lengths are Å. No random global state.
function may30_geometry(kind=:bcc;g=4)
    if kind===:bcc
        lattice=3.4Matrix{Float64}(I,3,3); basis=[0. .5;0. .5;0. .5]; reps=(g,g,g)
    elseif kind===:hcp
        a=3.2; lattice=[a a/2 0.;0 sqrt(3)*a/2 0.;0 0 sqrt(8/3)*a]
        basis=[0. 1/3;0. 1/3;0. .5]; reps=(5,5,3)
    elseif kind===:omega
        a=4.6; lattice=[a a/2 0.;0 sqrt(3)*a/2 0.;0 0 2.85]
        basis=[0. 1/3 2/3;0. 1/3 2/3;0. .5 .5]; reps=(4,4,5)
    else
        error("unknown fixture")
    end
    H=lattice*Diagonal(collect(reps))
    # Nonsymmetric homogeneous shear, shared by positions and cell.
    strain=[1.01 .023 -.011;0. .995 .017;0. 0 1.006]
    columns=[lattice*([i,j,k]+basis[:,b]) for b=axes(basis,2),i=0:reps[1]-1,j=0:reps[2]-1,k=0:reps[3]-1]
    X=reduce(hcat,vec(columns)); X=strain*X; H=strain*H
    for i in axes(X,2),a=1:3; X[a,i]+=.031sin(1.73i+.37a)+.009cos(.41i*a); end
    # These basis neighbors are checked explicitly in the comparison harness.
    z=[(22,13,23)[mod1(i,3)] for i in axes(X,2)]; z[1:2].=23
    # Independent exhaustive center-1 nearest image; assign V at both ends.
    vv_distance=Inf; vv_neighbor=0; vv_shift=(0,0,0)
    for j=2:length(z),sx=-1:1,sy=-1:1,sz=-1:1
        shift=(sx,sy,sz); distance=norm(X[:,j]-X[:,1]+H*collect(shift))
        if distance<vv_distance
            vv_distance=distance; vv_neighbor=j; vv_shift=shift
        end
    end
    z[vv_neighbor]=23
    (;positions=X,cell=H,species=z,vv_distance,vv_neighbor,vv_shift,kind=String(kind),seed="analytic sin/cos v1; no RNG")
end
