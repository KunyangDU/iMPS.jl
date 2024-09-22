function WeightProd(MPS::Vector,weight::Number)
    return (vcat([weight],ones(length(MPS)-1)) .* MPS)
end

function WeightSum(MPSs::Vector,weights::Vector,D_MPS::Int64)
    d = dims(domain(MPSs[1][1]))


    MPS = VariPlusMPS(WeightProd(MPSs[1],weights[1]), WeightProd(MPSs[2],weights[2]),d,D_MPS)

    for i in eachindex(MPSs)[3:end]
        MPS = VariPlusMPS(MPS,WeightProd(MPSs[i],weights[i]),d,D_MPS)
    end
    
    return MPS
end
