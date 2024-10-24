function addIntr!(Root::InteractionTreeNode,
    Opri::AbstractTensorMap,
    site::Int64,
    name::String,
    strength::Number,
    Z::Union{Nothing,AbstractTensorMap})

    addIntr1!(Root,Opri,site,name,strength,Z)
end

function addIntr!(Root::InteractionTreeNode,
    Opri::NTuple{1,AbstractTensorMap},
    site::NTuple{1,Int64},
    name::NTuple{1,String},
    strength::Number,
    Z::Union{Nothing,AbstractTensorMap})

    addIntr1!(Root,Opri[1],site[1],name[1],strength,Z)
end

function addIntr!(Root::InteractionTreeNode,
    Opri::NTuple{2,AbstractTensorMap},
    site::NTuple{2,Int64},
    name::NTuple{2,String},
    strength::Number,
    Z::Union{Nothing,AbstractTensorMap})

    addIntr2!(Root,Opri,site,name,strength,Z)
end

function addIntr!(Tree::InteractionTree,
    Opri::AbstractTensorMap,
    site::Int64,
    name::String,
    strength::Number,
    Z::Union{Nothing,AbstractTensorMap})

    addIntr1!(Tree.Root.children[1],Opri,site,name,strength,Z)
end

function addIntr!(Tree::InteractionTree,
    Opri::NTuple{1,AbstractTensorMap},
    site::NTuple{1,Int64},
    name::NTuple{1,String},
    strength::Number,
    Z::Union{Nothing,AbstractTensorMap})

    addIntr1!(Tree.Root.children[1],Opri[1],site[1],name[1],strength,Z)
end

function addIntr!(Tree::InteractionTree,
    Opri::NTuple{2,AbstractTensorMap},
    site::NTuple{2,Int64},
    name::NTuple{2,String},
    strength::Number,
    Z::Union{Nothing,AbstractTensorMap})

    addIntr2!(Tree.Root.children[1],Opri,site,name,strength,Z)
end

############# k ####################

function addIntr!(Root::InteractionTreeNode,
    Opri::AbstractTensorMap,
    Latt::AbstractLattice,k::Vector,
    name::String,
    strength::Number,
    string::Union{Nothing,AbstractTensorMap})
    L = size(Latt)
    for site in 1:L
        addIntr1!(Root,Opri,site,name,strength*exp(-1im*dot(k,coordinate(Latt,site))) / sqrt(L),string)
    end
end

function addIntr!(Tree::InteractionTree,
    Opri::AbstractTensorMap,
    Latt::AbstractLattice,k::Vector,
    name::String,
    strength::Number,
    string::Union{Nothing,AbstractTensorMap})
    addIntr!(Tree.Root.children[1],Opri,Latt,k,name,strength,string)
end
