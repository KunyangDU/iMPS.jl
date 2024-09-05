
function RightEnv(ψ::Vector,H::Vector,site::Int64=1)

    HR = RightMostEnv(ψ[L],H[L])

    for iL in L-1:-1:site+1
        HR = PushLeft(HR,ψ[iL],H[iL])
    end

    return HR
end

function RightMostEnv(ψi::AbstractTensorMap,Hi::AbstractTensorMap)
    @tensor HR[-1,-2,-3] ≔ ψi[-1,1]*Hi[1,2,-2]*ψi'[2,-3]
    HR = permute(HR,(1,),(2,3))
    return HR
end

function PushLeft(HR::AbstractTensorMap,ψi::AbstractTensorMap,Hi::AbstractTensorMap)
    @tensor tempHR[-1,-2,-3] ≔ ψi[-1,1,4]*Hi[1,5,2,-2]*ψi'[2,6,-3]*HR[4,5,6]
    HR = permute(tempHR,(1,),(2,3))
    return HR
end

function LeftEnv(ψ::Vector,H::Vector,site::Int64)
    HL = LeftMostEnv(ψ[1],H[1])

    for iL in 2:site-1
        HL = PushRight(HL,ψ[iL],H[iL])
    end

    return HL
end

function LeftMostEnv(ψi::AbstractTensorMap,Hi::AbstractTensorMap)

    @tensor HL[-1,-2,-3] ≔ ψi[-1,1]*Hi[1,-2,2]*ψi'[2,-3]
    HL = permute(HL,(1,2),(3,))

    return HL
end

function PushRight(HL::AbstractTensorMap,ψi::AbstractTensorMap,Hi::AbstractTensorMap)

    @tensor tempHL[-1,-2,-3] ≔ ψi[-1,4,1]*Hi[1,-2,2,5]*ψi'[6,2,-3]*HL[4,5,6]
    HL = permute(tempHL,(1,2),(3,))

    return HL
end