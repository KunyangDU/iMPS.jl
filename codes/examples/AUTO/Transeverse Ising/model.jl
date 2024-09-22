
function Hamiltonian(Latt::AbstractLattice;
    J::Number=1,h::Number=0,
    returntree::Bool=false)

    Root = InteractionTreeNode(IdentityOperator(0))

    LocalSpace = Spin2

    for i in 1:size(Latt)
        addIntr1!(Root,LocalSpace.Sz,"Sz",i,h)
    end
    
    for pair in neighbor(Latt)
        addIntr2!(Root,LocalSpace.SxSx,pair,("Sx","Sx"),J,nothing)
    end

    if returntree
        return InteractionTree(Root)
    else
        return AutomataMPO(InteractionTree(Root))
    end
    
end

function IsingRandMPS(Latt::AbstractLattice)
    return RandMPS(size(Latt);d=2)   
end

