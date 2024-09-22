
function Hamiltonian(Latt::AbstractLattice;
    Jx::Number=1,Jy::Number=Jx,Jz::Number=Jx,h::Number=0,
    returntree::Bool=false)

    Root = InteractionTreeNode(IdentityOperator(0))

    LocalSpace = Spin2

    for i in 1:size(Latt)
        addIntr1!(Root,LocalSpace.Sz,"Sz",i,h)
    end
    
    for pair in neighbor(Latt)
        addIntr2!(Root,LocalSpace.SxSx,pair,("Sx","Sx"),Jx,nothing)
        addIntr2!(Root,LocalSpace.SySy,pair,("Sy","Sy"),Jy,nothing)
        addIntr2!(Root,LocalSpace.SzSz,pair,("Sz","Sz"),Jz,nothing)
    end

    if returntree
        return InteractionTree(Root)
    else
        return AutomataMPO(InteractionTree(Root))  
    end
    
end

function HeisenbergRandMPS(Latt::AbstractLattice)
    return RandMPS(size(Latt);d=2)   
end

