

function EffHam(ψ::Vector,H::Vector,site::Int64)

    if site == 1
        HR = RightEnv(ψ,H,site)
        reduH = EffHam(H[site],HR)
    elseif site == length(ψ)
        HL = LeftEnv(ψ,H,site)
        reduH = EffHam(H[site],HL)
    else
        HR = RightEnv(ψ,H,site)
        HL = LeftEnv(ψ,H,site)
        reduH = EffHam(H[site],HL,HR)
    end

    return reduH
end

function EffHam(Hi::AbstractTensorMap{ComplexSpace,2,1},HR::AbstractTensorMap{ComplexSpace,1,2})
    @tensor reduH[-1,-2,-3,-4] ≔ Hi[-1,1,-3]*HR[-2,1,-4]
    reduH = permute(reduH,(1,2),(3,4))
    return reduH
end

function EffHam(Hi::AbstractTensorMap{ComplexSpace,1,2},HL::AbstractTensorMap{ComplexSpace,2,1})
    @tensor reduH[-1,-2,-3,-4] ≔ HL[-1,1,-3]*Hi[-2,-4,1]
    reduH = permute(reduH,(1,2),(3,4))
    return reduH
end

function EffHam(Hi::AbstractTensorMap{ComplexSpace,2,2},HL::AbstractTensorMap{ComplexSpace,2,1},HR::AbstractTensorMap{ComplexSpace,1,2})
    @tensor reduH[-1,-2,-3,-4,-5,-6] ≔ HL[-1,1,-4]*Hi[-2,2,-5,1]*HR[-3,2,-6]
    reduH = permute(reduH,(1,2,3),(4,5,6))
    return reduH
end


function EvolveOpr(Ham::AbstractTensorMap,τ::Number)
    return exp(-1im * Ham * τ)
end


function Apply(ψ::AbstractTensorMap{ComplexSpace,0,2},Opr::AbstractTensorMap{ComplexSpace,2,2})
    @tensor ψτ[-1,-2] ≔ ψ[1,2]*Opr[1,2,-1,-2]
    return permute(ψτ,(),(1,2))
end

function Apply(ψ::AbstractTensorMap{ComplexSpace,0,3},Opr::AbstractTensorMap{ComplexSpace,3,3})
    @tensor ψτ[-1,-2,-3] ≔ ψ[1,2,3]*Opr[1,2,3,-1,-2,-3]
    return permute(ψτ,(),(1,2,3))
end

function LeftApply(ψ::AbstractTensorMap{ComplexSpace,0,3},Opr::AbstractTensorMap{ComplexSpace,3,3})
    @tensor ψτ[-3,-1,-2] ≔ ψ[3,1,2]*Opr[1,2,3,-1,-2,-3]
    return permute(ψτ,(),(1,2,3))
end

#= function Apply(ψ::AbstractTensorMap{ComplexSpace,1,2},Opr::AbstractTensorMap{ComplexSpace,2,4})
    @tensor ψτ[-1,-2,-3] ≔ ψ[1,2,3]*Opr[1,-1,2,3,-2,-3]
    @show 127
    return permute(ψτ,(1,),(2,3))
end

function Apply(ψ::AbstractTensorMap{ComplexSpace,1,1},Opr::AbstractTensorMap{ComplexSpace,2,2})
    @tensor ψτ[-1,-2] ≔ ψ[1,2]*Opr[1,-1,2,-2]
    @show 133
    return permute(ψτ,(1,),(2,))
end

function LeftApply(ψ::AbstractTensorMap{ComplexSpace,0,2},Opr::AbstractTensorMap{ComplexSpace,2,2})
    @tensor ψτ[-1,-2] ≔ ψ[1,2]*Opr[1,2,-1,-2]
    @show 145
    return permute(ψτ,(),(1,2))
end

function EffHam0(EnvL::AbstractTensorMap{ComplexSpace,2,1},EnvR::AbstractTensorMap{ComplexSpace,1,2})
    @tensor reduH0[-1,-2,-3,-4] ≔ EnvL[-1,1,-3] * EnvR[-2,1,-4]
    return permute(reduH0,(1,2),(3,4))
end =#
