

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
            tempDict[subRoot.Opr.name] = QuantUniv(ψ, compress(canonicalize(AutomataMPO(InteractionTree(subRoot),L)))[1])
        end
        ObsDict[Root.Opr.name] = tempDict
    end
    Obsf = nothing
    return ObsDict
end


function calObs(Tree::ObserableTree,L::Int64 = treeheight(Tree.Root) - 1)
    Roots = Tree.Root.children
    ObsDict = Dict()
    for Root in Roots
        ObsDict[Root.Opr.name] = Dict()
    end
    return ObsDict
end

function calObs(
    ρ::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}},
    Obsf::ObserableForest,
    L::Int64=treeheight(Obsf.Roots) - 2;
    Z::Union{Nothing,Number} = nothing)
    Roots = Obsf.Roots.children
    ObsDict = Dict{String,Dict}()
    for Root in Roots
        tempDict = Dict{Tuple,Number}()
        for subRoot in Root.children
            cutparent!(subRoot)
            if isnothing(Z)
                tempDict[subRoot.Opr.name] = Trace(ρ, compress(canonicalize(AutomataMPO(InteractionTree(subRoot),L)))[1];Z=Trace(ρ))
            else
                tempDict[subRoot.Opr.name] = Trace(ρ, compress(canonicalize(AutomataMPO(InteractionTree(subRoot),L)))[1];Z=Z)
            end
        end
        ObsDict[Root.Opr.name] = tempDict
    end
    Obsf = nothing
    return ObsDict
end


