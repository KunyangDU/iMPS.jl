LocalSpace = Spin2

function Hamiltonian(Latt::AbstractLattice;
    J::Number=1,h::Number=0,
    returntree::Bool=false)

    Root = InteractionTreeNode()


    for i in 1:size(Latt)
        addIntr!(Root,LocalSpace.Sx,i,"Sx",h,nothing)
    end
    
    for pair in neighbor(Latt)
        addIntr!(Root,LocalSpace.SzSz,pair,("Sz","Sz"),J,nothing)
    end

    if returntree
        return InteractionTree(Root)
    else
        return AutomataMPO(InteractionTree(Root),size(Latt))  
    end
    
end

function InitialRandΨ(Latt::AbstractLattice)
    return RandMPS(size(Latt);d=LocalSpace.d)   
end

