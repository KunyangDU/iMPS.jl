

function reduHam(ψ::Vector,H::Vector,site::Int64)

    # 缩并顺序可以优化

    if site == 1
        HR = RightEnv(ψ,H,site)
        reduH = LeftMostReduHam(H[site],HR)
    elseif site == length(ψ)
        HL = LeftEnv(ψ,H,site)
        reduH = RightMostReduHam(H[site],HL)
    else
        HR = RightEnv(ψ,H,site)
        HL = LeftEnv(ψ,H,site)
        reduH = InnerReduHam(H[site],HL,HR)
    end

    return reduH
end

function LeftMostReduHam(Hi::AbstractTensorMap,HR::AbstractTensorMap)
    @tensor reduH[-1,-2,-3,-4] ≔ Hi[-1,1,-3]*HR[-2,1,-4]
    reduH = permute(reduH,(1,2),(3,4))
    return reduH
end

function RightMostReduHam(Hi::AbstractTensorMap,HL::AbstractTensorMap)
    @tensor reduH[-1,-2,-3,-4] ≔ HL[-1,1,-3]*Hi[-2,-4,1]
    reduH = permute(reduH,(1,2),(3,4))
    return reduH
end

function InnerReduHam(Hi::AbstractTensorMap,HL::AbstractTensorMap,HR::AbstractTensorMap)
    @tensor reduH[-1,-2,-3,-4,-5,-6] ≔ HL[-1,1,-4]*Hi[-2,2,-5,1]*HR[-3,2,-6]
    reduH = permute(reduH,(1,2,3),(4,5,6))
    return reduH
end
