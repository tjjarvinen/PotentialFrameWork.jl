function element_channel(model,sp)
    # Explicit element-only policy: reject isotopes and custom names. `==` in
    # AtomsBase uses wildcard matching, so strict identity is essential here.
    s=ChemicalSpecies(sp)
    z=Int(AtomsBase.atomic_number(s))
    isequal(s,ChemicalSpecies(z)) || throw(ArgumentError("ACE element channels do not accept isotopes or named species: $s"))
    i=findfirst(==(z),model.atomic_numbers)
    isnothing(i) && throw(ArgumentError("ACE model has no element channel for $s"))
    i
end

function prepare_environment!(env::FullEnvironment,sys,model::ResolvedModel;kwargs...)
    start=time_ns(); empty!(env.timings)
    channels=[element_channel(model,AtomsBase.species(sys,i)) for i=1:length(sys)]
    env.timings[:element_channels]=(time_ns()-start)*1e-9
    NL.prepare_full_environment!(env,sys,channels,length(model.atomic_numbers),NI.search_cutoff(model);kwargs...)
end

"""Owning AtomsBase extraction in numeric Å; H has cell vectors in columns.
Open lattice generators are zeroed before numerical evaluation.
"""
function ace_inputs(sys,model;environment=FullEnvironment(),kwargs...)
    start=time_ns()
    n=length(sys); cv=AtomsBase.cell_vectors(sys); pb=AtomsBase.periodicity(sys)
    X=reshape([ustrip(u"Å",AtomsBase.position(sys,i)[a]) for a=1:3,i=1:n],3,n)
    H=[pb[j] ? ustrip(u"Å",cv[j][a]) : zero(ustrip(u"Å",cv[j][a])) for a=1:3,j=1:3]
    elapsed=(time_ns()-start)*1e-9
    prepare_environment!(environment,sys,model;kwargs...)
    environment.timings[:boundary_extraction]=elapsed
    (;positions=X,cell=H,environment)
end
