function addObs!(Tree::ObserableTree{N},
    Opri::AbstractTensorMap,
    site::Int64,
    name::String,
    Z::Union{Nothing,AbstractTensorMap};ObsName = nothing
    ) where N 
    isnothing(ObsName) && (ObsName = name)
    ind = findfirst(map(x -> x.Opr.name == ObsName,Tree.Root.children))
    if isnothing(ind)
        _addBranch!(Obs,ObsName)
        addIntr!(Tree.Root.children[end],Opri,site,name,1,Z)
    else
        addIntr!(Tree.Root.children[ind],Opri,site,name,1,Z)
    end
end

function addObs!(Tree::ObserableTree{N},
    Opri::NTuple{2,AbstractTensorMap},
    site::NTuple{2,Int64},
    name::NTuple{2,String},
    Z::Union{Nothing,AbstractTensorMap};ObsName = nothing
    ) where N 
    isnothing(ObsName) && (ObsName = string(name...))
    ind = findfirst(map(x -> x.Opr.name == ObsName,Tree.Root.children))
    if isnothing(ind)
        _addBranch!(Obs,ObsName)
        addIntr!(Tree.Root.children[end],Opri,site,name,1,Z)
    else
        addIntr!(Tree.Root.children[ind],Opri,site,name,1,Z)
    end
end

function addObs!(Obsf::ObserableForest{N},
    Opri::Union{AbstractTensorMap,NTuple{2,AbstractTensorMap}},
    site::Union{Int64,NTuple{2,Int64}},
    name::Union{String,NTuple{2,String}},
    Z::Union{Nothing,AbstractTensorMap};ObsName = nothing
    ) where N
    isnothing(ObsName) && (ObsName = string(name...))
    ind = findfirst(map(x -> x.Opr.name == ObsName,Obsf.Roots.children))
    if isnothing(ind)
        _addBranch!(Obsf,ObsName)
        _addBranch!(Obsf.Roots.children[end],Tuple(site))
        addIntr!(Obsf.Roots.children[end].children[end],Opri,site,name,1,Z)
    else
        indt = findfirst(map(x -> x.Opr.name == Tuple(site),Obsf.Roots.children[ind].children))
        if isnothing(indt)
            _addBranch!(Obsf.Roots.children[end],Tuple(collect(site)))
            addIntr!(Obsf.Roots.children[end].children[end],Opri,site,name,1,Z)
        else
            @error "site already exist"
        end
    end
end

function _addBranch!(Obs::ObserableTree,name::String)
    addchild!(Obs.Root,InteractionTreeNode(IdentityOperator(0,name)))
end

function _addBranch!(Obsf::ObserableForest,name::Union{String,Tuple})
    addchild!(Obsf.Roots,InteractionTreeNode(IdentityOperator(0,name)))
end

function _addBranch!(Root::InteractionTreeNode,name::Union{String,Tuple})
    addchild!(Root,InteractionTreeNode(IdentityOperator(0,name)))
end


