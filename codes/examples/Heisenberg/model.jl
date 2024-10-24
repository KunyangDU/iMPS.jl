LocalSpace = Spin2

function Hamiltonian(Latt::AbstractLattice;
    Jx::Number=1,Jy::Number=Jx,Jz::Number=Jx,h::Number=0,
    returntree::Bool=false)

    Root = InteractionTreeNode()

<<<<<<< HEAD
    I = diagm(ones(2))
    I0 = zeros(2,2)
    σx = [0 1;1 0]
    σy = [0 -1im;1im 0]
    σz = [1 0;0 -1]

    d = 2
    D = 5

    MPO = Vector{AbstractTensorMap{ComplexSpace,2,2}}(undef, L)

    idt = ℂ^1
    phys = (ℂ^d)'
    bond = ℂ^D
    for i in 1:L
        if i == 1
            data = reshape([-hz*σz -Jx*σx -Jy*σy -Jz*σz I],d,D,d,1)
            M = BlockMPO(data,phys,idt,phys,bond)
        elseif i == L
            data = reshape([I;σx;σy;σz;-hz*σz],d,D,d,1)
            M = BlockMPO(data,phys,bond,phys,idt)
        else
            data = reshape([
                I I0 I0 I0 I0
                σx I0 I0 I0 I0
                σy I0 I0 I0 I0
                σz I0 I0 I0 I0
                -hz*σz -Jx*σx -Jy*σy -Jz*σz I
            ],d,D,d,D)
            M = BlockMPO(data,phys,bond,phys,bond)
        end
        MPO[i] = M
=======
    for i in 1:size(Latt)
        addIntr!(Root,LocalSpace.Sz,i,"Sz",h,nothing)
>>>>>>> 8ea8417fd317c4adb4f58a9cd6b4c299e7c2f40e
    end
    
    for pair in neighbor(Latt)
        addIntr!(Root,LocalSpace.SxSx,pair,("Sx","Sx"),Jx,nothing)
        addIntr!(Root,LocalSpace.SySy,pair,("Sy","Sy"),Jy,nothing)
        addIntr!(Root,LocalSpace.SzSz,pair,("Sz","Sz"),Jz,nothing)
    end

<<<<<<< HEAD
    MPO = Vector{AbstractTensorMap{ComplexSpace,2,2}}(undef,size(Latt))
=======
    if returntree
        return InteractionTree(Root)
    else
        return AutomataMPO(InteractionTree(Root),size(Latt))  
    end
>>>>>>> 8ea8417fd317c4adb4f58a9cd6b4c299e7c2f40e
    
end

function InitialRandΨ(Latt::AbstractLattice)
    return RandMPS(size(Latt);d=LocalSpace.d)   
end


function SkMPOs(Latt::AbstractLattice,k::Vector)
    MPOs = map(x -> compress(canonicalize(KOprMPO(x[1],Latt,k,x[2],nothing))),[(LocalSpace.Sx,"Sxk"),(LocalSpace.Sy,"Syk"),(LocalSpace.Sz,"Szk")])
    return map(x -> [mpo[x] for mpo in MPOs],1:2)
end
