
LocalSpace = SpinlessFermion
function Hamiltonian(Latt::AbstractLattice;t::Number=1,U::Number=0,μ::Number=0,returntree::Bool=false)

    Root = InteractionTreeNode(IdentityOperator(0))

    for i in 1:size(Latt)
        addIntr!(Root,LocalSpace.n,i,"n",-μ,nothing)
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

function InitialRandΨ(Latt::AbstractLattice)
    return RandMPS(size(Latt);d=LocalSpace.d)   
end

function ParticleNumber(Latt::AbstractLattice,site::Int64)
    Root = InteractionTreeNode(IdentityOperator(0))

    addIntr!(Root,LocalSpace.n,site,"n",1,nothing)

    return AutomataMPO(InteractionTree(Root))
end

function ParticleNumber(Latt::AbstractLattice)
    Root = InteractionTreeNode(IdentityOperator(0))

    for site in 1:size(Latt)
        addIntr!(Root,LocalSpace.n,site,"n",1,nothing)
    end

    return AutomataMPO(InteractionTree(Root)) 
end


