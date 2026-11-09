# Deterministic perturbed beta-Ti/B2-like bcc cells; construction is never timed.
function bcc_geometry(g,z;spacing=3.4)
    n=2g^3
    X=[begin
        cell=(i-1)÷2; shift=isodd(i) ? 0. : .5
        spacing*((a==1 ? cell%g : a==2 ? (cell÷g)%g : cell÷g^2)+shift)+.025*sin(i*(a+.31))
    end for a=1:3,i=1:n]
    H=[a==b ? g*spacing : 0. for a=1:3,b=1:3]
    numbers=Int.(z)
    alloy=length(z)==3 && Set(numbers)==Set([13,22,23])
    labels=if alloy
        # Approximately Ti-6Al-4V in atomic fractions, distributed deterministically.
        [z[findfirst(==((37*(i-1))%100<10 ? 13 : (37*(i-1))%100<14 ? 23 : 22),numbers)] for i=1:n]
    else
        [z[mod1(i,length(z))] for i=1:n]
    end
    (;positions=X,cell=H,species=labels)
end
