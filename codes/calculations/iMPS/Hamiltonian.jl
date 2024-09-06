

function ReduHam1(ψ::Vector,H::Vector,site::Int64)

    if site == 1
        HR = RightEnv(ψ,H,site)
        reduH = ReduHam1(H[site],HR)
    elseif site == length(ψ)
        HL = LeftEnv(ψ,H,site)
        reduH = ReduHam1(H[site],HL)
    else
        HR = RightEnv(ψ,H,site)
        HL = LeftEnv(ψ,H,site)
        reduH = ReduHam1(H[site],HL,HR)
    end

    return reduH
end

function ReduHam1(Hi::AbstractTensorMap{ComplexSpace,2,1},HR::AbstractTensorMap{ComplexSpace,1,2})
    @tensor reduH[-1,-2,-3,-4] ≔ Hi[-1,1,-3]*HR[-2,1,-4]
    reduH = permute(reduH,(1,2),(3,4))
    return reduH
end

function ReduHam1(Hi::AbstractTensorMap{ComplexSpace,1,2},HL::AbstractTensorMap{ComplexSpace,2,1})
    @tensor reduH[-1,-2,-3,-4] ≔ HL[-1,1,-3]*Hi[-2,-4,1]
    reduH = permute(reduH,(1,2),(3,4))
    return reduH
end

function ReduHam1(Hi::AbstractTensorMap{ComplexSpace,2,2},HL::AbstractTensorMap{ComplexSpace,2,1},HR::AbstractTensorMap{ComplexSpace,1,2})
    @tensor reduH[-1,-2,-3,-4,-5,-6] ≔ HL[-1,1,-4]*Hi[-2,2,-5,1]*HR[-3,2,-6]
    reduH = permute(reduH,(1,2,3),(4,5,6))
    return reduH
end
