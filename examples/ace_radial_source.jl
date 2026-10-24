# Self-contained trusted source bundle; no captured globals or package installs.
struct PortableRadial end
function (::PortableRadial)(v,d,r,ps,config)
    rc=config.cutoff; x=r/rc; q=max(one(r)-x,zero(r))
    q2=q*q; q3=q2*q; q4=q2*q2
    rate=ps.rate[1]; amplitude=ps.amplitude[1]; ex=exp(-rate*r)
    for i in eachindex(v)
        n=config.powers[i]
        xn=one(r); dx=zero(r)
        for j=1:n
            dx=dx*x+xn/rc
            xn*=x
        end
        v[i]=amplitude*xn*q4*ex
        isnothing(d) || (d[i]=amplitude*ex*(dx*q4-4xn*q3/rc-rate*xn*q4))
    end
    nothing
end
radial_factory(config,metadata)=PortableRadial()
