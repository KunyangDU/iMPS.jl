# Contract
# Pure merge between two local MPO or local MPS

function Contract(
    Opr1::AbstractTensorMap{ComplexSpace,1,3},
    Opr2::AbstractTensorMap{ComplexSpace,2,2})
    @tensor tempOpr[-1,-2,-3,-4,-5,-6] ≔ Opr1[-2,-3,-4,1]*Opr2[-1,1,-5,-6]
    return permute(tempOpr,(1,2),(3,4,5,6))
end

function Contract(
    Opr1::AbstractTensorMap{ComplexSpace,2,2},
    Opr2::AbstractTensorMap{ComplexSpace,1,3})
    @tensor tempOpr[-1,-2,-3,-4,-5,-6] ≔ Opr1[1,-2,-3,-4]*Opr2[-1,1,-5,-6]
    return permute(tempOpr,(1,2),(3,4,5,6))
end

function Contract(ψ1::AbstractTensorMap{ComplexSpace,1,2},
    ψ2::AbstractTensorMap{ComplexSpace,0,3})

    @tensor ψm[-1,-2,-3,-4] ≔ ψ1[1,-1,-2] * ψ2[1,-3,-4]
    return permute(ψm,(),(1,2,3,4))
end

function Contract(ψ1::AbstractTensorMap{ComplexSpace,0,3},
    ψ2::AbstractTensorMap{ComplexSpace,1,2})

    @tensor ψm[-1,-2,-3,-4] ≔ ψ1[-1,-2,1] * ψ2[1,-3,-4]
    return permute(ψm,(),(1,2,3,4))
end

function InnerProd(ψ₁::Vector,ψ₂::Vector)
    EnvL = InitialLeftEnv(;order=2)
    EnvR = RightEnv(ψ₁,ψ₂,1)
    innerprod = @tensor EnvL[1,3]*ψ₁[1][1,5,2]*ψ₂[1]'[3,5,4]*EnvR[2,4]
    return innerprod[1]
end

function Apply(ψi::AbstractTensorMap{ComplexSpace,0,2},Opr::AbstractTensorMap{ComplexSpace,2,2})
    @tensor ψia[-1,-2] ≔ ψi[1,2]*Opr[1,2,-1,-2]
    return permute(ψia,(),(1,2))
end

function Apply(ψi::AbstractTensorMap{ComplexSpace,0,3},Opr::AbstractTensorMap{ComplexSpace,3,3})
    @tensor ψia[-1,-2,-3] ≔ ψi[1,3,2]*Opr[1,3,2,-1,-2,-3]
    return permute(ψia,(),(1,2,3))
end

function Apply(ψi::AbstractTensorMap{ComplexSpace,0,4},Opr::AbstractTensorMap{ComplexSpace,4,4})
    @tensor ψia[-1,-2,-3,-4] ≔ ψi[1,3,4,2]*Opr[1,3,4,2,-1,-2,-3,-4]
    return permute(ψia,(),(1,2,3,4))
end

