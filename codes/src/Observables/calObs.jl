

function calObs(ψ::Vector{Union{AbstractTensorMap{ComplexSpace,1,2},AbstractTensorMap{ComplexSpace,0,3}}},
    Opr::Vector{AbstractTensorMap{ComplexSpace,2,2}})
    return QuantUniv(ψ,Opr)
end

function calObs(Tree::ObserableTree,L::Int64 = treeheight(Tree.Root) - 1)
    Roots = Tree.Root.children
    ObsDict = Dict()
    for Root in Roots
        ObsDict[Root.Opr.name] = Dict()
    end
    return ObsDict
end

function calObs(ψ::Vector{Union{AbstractTensorMap{ComplexSpace,1,2},AbstractTensorMap{ComplexSpace,0,3}}},
    Obsf::ObserableForest,L::Int64=treeheight(Obsf.Roots) - 2)
    Roots = Obsf.Roots.children
    ObsDict = Dict{String,Dict}()
    for Root in Roots
        tempDict = Dict{Tuple,Number}()
        for subRoot in Root.children
            cutparent!(subRoot)
            tempDict[subRoot.Opr.name] = QuantUniv(ψ, AutomataMPO(InteractionTree(subRoot),L))
        end
        ObsDict[Root.Opr.name] = tempDict
    end
    Obsf = nothing
    return ObsDict
end


