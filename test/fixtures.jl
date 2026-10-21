using AtomsBase, Unitful
# A separate AtomsBase implementation ensures no dependency on a concrete helper system.
struct TestSystem{T} <: AtomsBase.AbstractSystem{3}
    xyz::Matrix{T}
    H::Matrix{T}
    pbc::NTuple{3,Bool}
    labels::Vector{ChemicalSpecies}
    TestSystem(xyz::Matrix{T}, H::Matrix{T}, pbc, labels) where T =
        new{T}(xyz,H,pbc,ChemicalSpecies.(labels))
end
Base.length(s::TestSystem)=size(s.xyz,2)
AtomsBase.position(s::TestSystem,i::Integer)=s.xyz[:,i]*u"Å"
AtomsBase.cell_vectors(s::TestSystem)=ntuple(a -> s.H[:,a]*u"Å",3)
AtomsBase.periodicity(s::TestSystem)=s.pbc
AtomsBase.species(s::TestSystem,i::Integer)=s.labels[i]

function fixture(n; L=9., pbc=(true,true,true), cluster=false, mixed=true)
    xyz=[L*mod(sin(i*(a+0.71))*437.13,1) for a=1:3,i=1:n]
    cluster && (xyz[:,1:cld(n,2)] .*= 0.18)
    H=zeros(3,3)
    for a=1:3; H[a,a]=L; end
    labels=[mixed && isodd(i) ? :He : :H for i=1:n]
    TestSystem(xyz,H,pbc,labels)
end
potential()=SmoothPair([:H,:He],[1.0 0.7; 0.7 1.3],2.0)

# Independent O(N²) oracle in original order, no production traversal or masks.
function oracle(s,p)
    n=length(s); f=zeros(3,n); v=zeros(3,3); energy=0.0; pairs=Set{Tuple{Int,Int}}()
    for i=1:n, j=i+1:n
        d=s.xyz[:,j]-s.xyz[:,i]
        for a=1:3
            if s.pbc[a]
                L=s.H[a,a]
                # deliberately different implementation of minimum-image mapping
                while d[a]>L/2; d[a]-=L; end
                while d[a]<-L/2; d[a]+=L; end
            end
        end
        r2=sum(abs2,d)
        r2<p.radius^2 || continue
        push!(pairs,(i,j))
        a=p.amplitudes[findfirst(isequal(s.labels[i]),p.labels),findfirst(isequal(s.labels[j]),p.labels)]
        q=1-r2/p.radius^2
        energy+=a*q^3
        force=(-6*a*q^2/p.radius^2)*d
        f[:,i]+=force; f[:,j]-=force
        v-=d*transpose(force)
    end
    (;energy,forces=f,virial=v,pairs)
end
