
function SETTN(
    H::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}},
    β::Number,order::Int64,D_MPO::Int64)
    d = dim(codomain(H[1]))
    L = length(H)

    tempρ = IdentityMPO(L,d)

    ρ = deepcopy(tempρ)

    for i in 1:order 
        tempρ = VariProdMPO(H,tempρ,d,D_MPO)
        ρ = let 
            temp = deepcopy(tempρ)
            temp[1] = SETTNCoeff(β,i) * temp[1]
            VariPlusMPO(ρ,temp,d,D_MPO)
        end
    end
    ρ[1] = ρ[1]/norm(ρ[1])

    return ρ
end

function SETTNCoeff(β::Number,order::Int64)
    return (-β)^order / factorial(order)
end