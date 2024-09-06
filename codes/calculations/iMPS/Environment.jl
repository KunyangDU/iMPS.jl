#修改
function RightEnv(ψ::Vector,H::Vector,site::Int64=1)

    HR = RightEnv(ψ[L],H[L])

    for iL in L-1:-1:site+1
        HR = PushLeft(HR,ψ[iL],H[iL])
    end

    return HR
end

function RightEnv(ψi::AbstractTensorMap{ComplexSpace,1,1},Hi::AbstractTensorMap{ComplexSpace,1,2})
    @tensor HR[-1,-2,-3] ≔ ψi[-1,1]*Hi[1,2,-2]*ψi'[2,-3]
    return permute(HR,(1,),(2,3))
end

function PushLeft(HR::AbstractTensorMap{ComplexSpace,1,2},ψi::AbstractTensorMap{ComplexSpace,1,2},Hi::AbstractTensorMap{ComplexSpace,2,2})
    @tensor tempHR[-1,-2,-3] ≔ ψi[-1,1,4]*Hi[1,5,2,-2]*ψi'[2,6,-3]*HR[4,5,6]
    return permute(tempHR,(1,),(2,3))
end

function LeftEnv(ψ::Vector,H::Vector,site::Int64)
    HL = LeftEnv(ψ[1],H[1])

    for iL in 2:site-1
        HL = PushRight(HL,ψ[iL],H[iL])
    end

    return HL
end

function LeftEnv(ψi::AbstractTensorMap{ComplexSpace,1,1},Hi::AbstractTensorMap{ComplexSpace,2,1})
    @tensor HL[-1,-2,-3] ≔ ψi[-1,1]*Hi[1,-2,2]*ψi'[2,-3]
    return permute(HL,(1,2),(3,))
end

function PushRight(HL::AbstractTensorMap{ComplexSpace,2,1},ψi::AbstractTensorMap{ComplexSpace,1,2},Hi::AbstractTensorMap{ComplexSpace,2,2})
    @tensor tempHL[-1,-2,-3] ≔ ψi[-1,4,1]*Hi[1,-2,2,5]*ψi'[6,2,-3]*HL[4,5,6]
    return permute(tempHL,(1,2),(3,))
end