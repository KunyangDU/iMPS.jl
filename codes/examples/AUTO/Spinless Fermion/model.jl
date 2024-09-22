

function Hamiltonian(Latt::AbstractLattice;t::Number=1,U::Number=0,μ::Number=0,returntree::Bool=false)

    Root = InteractionTreeNode(IdentityOperator(0))

    LocalSpace = SpinlessFermion

    for i in 1:size(Latt)
        addIntr1!(Root,SpinlessFermion.n,"n",i,-μ)
    end
    
    for pair in neighbor(Latt)
        addIntr2!(Root,LocalSpace.FFdag,pair,("F","F⁺"),t,LocalSpace.Z)
        addIntr2!(Root,LocalSpace.FdagF,pair,("F⁺","F"),t,LocalSpace.Z)
    end

    if returntree
        return InteractionTree(Root)
    else
        return AutomataMPO(InteractionTree(Root))  
    end
end

function SpinlessFermionRandMPS(Latt::AbstractLattice)
    return RandMPS(size(Latt);d=2)   
end

function ParticleNumber(Latt::AbstractLattice,site::Int64)
    Root = InteractionTreeNode(IdentityOperator(0))

    LocalSpace = SpinlessFermion

    addIntr1!(Root,LocalSpace.n,"n",site,1)

    return AutomataMPO(InteractionTree(Root))
end

function ParticleNumber(Latt::AbstractLattice)
    Root = InteractionTreeNode(IdentityOperator(0))

    LocalSpace = SpinlessFermion

    for i in 1:size(Latt)
        addIntr1!(Root,LocalSpace.n,"n",i,1)
    end

    return AutomataMPO(InteractionTree(Root)) 
end


