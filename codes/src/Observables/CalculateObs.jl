

function calObs(ψ::Vector{Union{AbstractTensorMap{ComplexSpace,1,2},AbstractTensorMap{ComplexSpace,0,3}}},
    Opr::Vector{AbstractTensorMap{ComplexSpace,2,2}})
    return QuantUniv(ψ,Opr)
end



