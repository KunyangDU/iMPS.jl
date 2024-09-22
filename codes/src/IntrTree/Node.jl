using AbstractTrees

mutable struct InteractionTreeNode

    Opr::Union{Nothing,AbstractLocalOperator}
    parent::Union{Nothing,InteractionTreeNode}
    children::Vector{InteractionTreeNode}

    function InteractionTreeNode(
        Opr::Union{Nothing,AbstractLocalOperator},
        parent::InteractionTreeNode,
        children::Vector{InteractionTreeNode}=InteractionTreeNode[],
    )
    return new(Opr,parent,children)
    end

    function InteractionTreeNode(
        Opr::Union{Nothing,AbstractLocalOperator},
        children::Vector{InteractionTreeNode}=InteractionTreeNode[],
    )
    return new(Opr,nothing,children)
    end

    InteractionTreeNode() = InteractionTreeNode("O", nothing)
end

AbstractTrees.nodevalue(node::InteractionTreeNode) = node.Opr
AbstractTrees.parent(node::InteractionTreeNode) = node.parent
AbstractTrees.children(node::InteractionTreeNode) = node.children
AbstractTrees.ParentLinks(::Type{InteractionTreeNode}) = StoredParents()
AbstractTrees.ChildIndexing(::Type{InteractionTreeNode}) = IndexedChildren()
AbstractTrees.NodeType(::Type{InteractionTreeNode}) = HasNodeType()
AbstractTrees.nodetype(::Type{InteractionTreeNode}) = InteractionTreeNode

function addchild!(node::InteractionTreeNode, child::InteractionTreeNode)
    isnothing(child.parent) ? child.parent = node : @assert child.parent == node
    push!(node.children, child)
    return nothing
end

function addchild!(node::InteractionTreeNode, Opr::AbstractLocalOperator)
    addchild!(node,InteractionTreeNode(Opr))
    return nothing
end


function Base.show(io::IO,Root::InteractionTreeNode)
    print_tree(Root;maxdepth = 16)
    return nothing
end


struct InteractionTree{N}
    Root::InteractionTreeNode
    function InteractionTree(child::InteractionTreeNode...)
         N = length(child)
         
         Root = InteractionTreeNode(nothing)
         for i in child
              addchild!(Root, i)
         end
         return new{N}(Root)
    end
end

