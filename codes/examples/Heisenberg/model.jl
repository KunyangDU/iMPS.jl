LocalSpace = Spin2

function Hamiltonian(Latt::AbstractLattice;
    Jx::Number=1,Jy::Number=Jx,Jz::Number=Jx,h::Number=0,
    returntree::Bool=false)

    Root = InteractionTreeNode()

    for i in 1:size(Latt)
        addIntr!(Root,LocalSpace.Sz,i,"Sz",h,nothing)
    end
    
    for pair in neighbor(Latt)
        addIntr!(Root,LocalSpace.SxSx,pair,("Sx","Sx"),Jx,nothing)
        addIntr!(Root,LocalSpace.SySy,pair,("Sy","Sy"),Jy,nothing)
        addIntr!(Root,LocalSpace.SzSz,pair,("Sz","Sz"),Jz,nothing)
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


function SkMPOs(Latt::AbstractLattice,k::Vector)
    MPOs = map(x -> compress(canonicalize(KOprMPO(x[1],Latt,k,x[2],nothing))),[(LocalSpace.Sx,"Sxk"),(LocalSpace.Sy,"Syk"),(LocalSpace.Sz,"Szk")])
    return map(x -> [mpo[x] for mpo in MPOs],1:2)
end
