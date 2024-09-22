function addIntr1!(Root::InteractionTreeNode,Opri::AbstractTensorMap,name::String,site::Int64,strength::Number)
    tempOpr = LocalOperator(Opri,name,site)
    addIntr1!(Root,tempOpr,strength)
end

function addIntr1!(Root::InteractionTreeNode,Opr::LocalOperator,strength::Number)
    current_node = Root
    current_site = 1

    # add the identity
    while current_site < Opr.site

        tempIdOpr = IdentityOperator(getIdTensor(Opr),current_site)

        indId = findfirst(x -> isequal(x.Opr,tempIdOpr),current_node.children)
        if isnothing(indId)
            addchild!(current_node,tempIdOpr)
            current_node = current_node.children[end]
        else
            current_node = current_node.children[indId]
        end

        current_site += 1
    end

    # add the onsite Opr 
    addchild!(current_node, Opr)
    current_node.children[end].Opr.strength = strength

end