
function AutomataMPO(Tree::InteractionTree,L::Int64 = treeheight(Tree.Root) - 1)
    Root = Tree.Root

    MPO = let 
        tempMPO = Vector{AbstractTensorMap{ComplexSpace,2,2}}(undef,L)

        last_leaves = []
        last_roots = Root.children

        idtensor = getIdTensor(last_roots[1].children[1].Opr)
        
        last_inverse_root = 0
        next_inverse_root = 0
        
        for iL in 1:L

            if next_inverse_root == 0 && !isempty(findall(x -> isempty(x.children),vcat([lastroot.children for lastroot in last_roots]...)))
                next_inverse_root = 1
            end

            next_leaves = []
            leaves_inds = []

            next_roots = []
            roots_inds = []

            for (lastind,last_root) in enumerate(last_roots)
                for next_subtree in last_root.children
                    if isempty(next_subtree.children)
                        push!(next_leaves,next_subtree)
                        push!(leaves_inds,(next_inverse_root,lastind + last_inverse_root,length(next_leaves)))
                    else
                        push!(next_roots,next_subtree)
                        push!(roots_inds,(length(next_roots) + next_inverse_root,lastind + last_inverse_root,length(next_roots)))
                    end
                end
            end

            localMPOdims = length.((next_roots,last_roots)) .+ (next_inverse_root,last_inverse_root)
            localMPO = FillBlockTensor((1,1),localMPOdims,last_inverse_root*idtensor)

            for inds in leaves_inds
                localMPO += FillBlockTensor(inds[1:2],localMPOdims,let 
                    localOpr = next_leaves[inds[3]].Opr.Opri
                    strength = next_leaves[inds[3]].Opr.strength
                    if isnan(strength)
                        localOpr
                    else
                        localOpr*strength
                    end
                end)
            end

            for inds in roots_inds
                localMPO += FillBlockTensor(inds[1:2],localMPOdims,let 
                    localOpr = next_roots[inds[3]].Opr.Opri
                    strength = next_roots[inds[3]].Opr.strength
                    if isnan(strength)
                        localOpr
                    else
                        localOpr*strength
                    end
                end)
            end
            
            last_leaves = next_leaves
            last_roots = next_roots
            last_inverse_root = next_inverse_root

            tempMPO[iL] = localMPO
        end

        tempMPO
    end
    
    return MPO
end

function FillBlockTensor(inds::NTuple{2,Int64},dims::NTuple{2,Int64},tblock::AbstractTensorMap)
    matr = zeros(dims...)
    matr[inds...] = 1
    return TensorMap(matr,map(x -> ℂ^x,dims)...) ⊗ tblock
end

