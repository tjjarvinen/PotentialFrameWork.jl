# Narrow legacy ACE1 orthogonal-polynomial family, represented as an ordinary
# explicit analytic callable. Mathematical recurrence/config only; no JuLIP runtime.
struct LegacyPolynomial end
function (::LegacyPolynomial)(v,d,r,ps,c)
    fill!(v,zero(r)); isnothing(d) || fill!(d,zero(r))
    r>=c.cutoff && return nothing
    s=(r-c.rin)/(c.r0-c.rin); h=c.q-c.p
    y=one(r); dy=zero(r)
    if r>c.rin
        den=one(r)+s^h; u=s^c.q/den
        du=(c.q*s^(c.q-1)/den-(h==0 ? zero(r) : h*s^(c.q+h-1)/den^2))/(c.r0-c.rin)
        y=inv(one(r)+c.a*u); dy=-c.a*du*y*y
    end
    x=c.y1+(y-c.x1)*(c.y2-c.y1)/(c.x2-c.x1)
    dx=dy*(c.y2-c.y1)/(c.x2-c.x1)
    ((c.pl>0 && x<c.tl) || (c.pr>0 && x>c.tr)) && return nothing
    left=x-c.tl; right=x-c.tr
    taper=left^c.pl*right^c.pr
    dtaper=(c.pl>0 ? c.pl*left^(c.pl-1)*right^c.pr : zero(r))+
           (c.pr>0 ? c.pr*left^c.pl*right^(c.pr-1) : zero(r))
    env=one(r); denv=zero(r)
    if c.envelope_power>0
        p=c.envelope_power; sr=r/c.envelope_r0; sc=c.cutoff/c.envelope_r0
        env=sr^(-p)-sc^(-p)+p*sc^(-p-1)*(sr-sc)
        denv=p/c.envelope_r0*(sc^(-p-1)-sr^(-p-1))
    end
    prev=zero(r); dprev=zero(r); value=c.A[1]*taper; deriv=c.A[1]*dtaper*dx
    for n in eachindex(c.A)
        if n>1
            t=c.A[n]*x+c.B[n]
            next=t*value+c.C[n]*prev
            dnext=t*deriv+c.A[n]*dx*value+c.C[n]*dprev
            prev=value; dprev=deriv; value=next; deriv=dnext
        end
        index=n+c.offset
        v[index]=value*env
        isnothing(d) || (d[index]=deriv*env+value*denv)
    end
    nothing
end
