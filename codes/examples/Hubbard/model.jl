LocalSpace = Spin2Fermion

function Hamiltonian(Latt::AbstractLattice;t::Number=1,U::Number=0,μ::Number=0,returntree::Bool=false)

    Root = InteractionTreeNode(IdentityOperator(0))

    for i in 1:size(Latt)
        addIntr!(Root,LocalSpace.n,i,"n",-μ,nothing)
        addIntr!(Root,LocalSpace.nd,i,"nd",U,nothing)
    end
    
    for pair in neighbor(Latt)
        addIntr!(Root,LocalSpace.FFdagUp,pair,("F↑","F⁺↑"),t,Spin2Fermion.Z)
        addIntr!(Root,LocalSpace.FFdagDown,pair,("F↓","F⁺↓"),t,Spin2Fermion.Z)
        addIntr!(Root,LocalSpace.FdagFUp,pair,("F⁺↑","F↑"),t,Spin2Fermion.Z)
        addIntr!(Root,LocalSpace.FdagFDown,pair,("F⁺↓","F↓"),t,Spin2Fermion.Z)
    end

    if returntree
        return InteractionTree(Root)
    else
        return AutomataMPO(InteractionTree(Root),size(Latt))  
    end

end

function ParticleNumber(Latt::AbstractLattice,site::Int64)
    Root = InteractionTreeNode(IdentityOperator(0))

    addIntr!(Root,LocalSpace.n,site,"n",1,nothing)

    return InteractionTree(Root)
end

function ParticleNumber(Latt::AbstractLattice)
    Root = InteractionTreeNode(IdentityOperator(0))

    LocalSpace = Spin2Fermion

    for i in 1:size(Latt)
        addIntr!(Root,LocalSpace.n,i,"n",1,nothing)
    end

    return InteractionTree(Root)   
end

function cupMPO(Latt::AbstractLattice,k::Vector)
    MPOs = map(y -> [mpo[y] for mpo in let 
    map(x -> compress(canonicalize(KOprMPO(x[1],Latt,k,x[2],LocalSpace.Z))),[(LocalSpace.Fup,"Fupk"),(LocalSpace.Fdagup,"Fupdagk")])
    end ],1:2)
    return MPOs
end

function cdownMPO(Latt::AbstractLattice,k::Vector)
    MPOs = map(y -> [mpo[y] for mpo in let 
    map(x -> compress(canonicalize(KOprMPO(x[1],Latt,k,x[2],LocalSpace.Z))),[(LocalSpace.Fdown,"Fdownk"),(LocalSpace.Fdagdown,"Fdowndagk")])
    end ],1:2)
    return MPOs
end

function cMPO(Latt::AbstractLattice,k::Vector)
    cupMPOs = cupMPO(Latt,k)
    cdownMPOs = cdownMPO(Latt,k)
    return map(x -> [cupMPOs[x],cdownMPOs[x]],1:2)
end


