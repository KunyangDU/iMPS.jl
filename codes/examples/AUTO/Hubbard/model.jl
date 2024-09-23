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


