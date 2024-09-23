
# MERGE
# along with some special operations like SVD or Environment

function RightMerge(Σ::AbstractTensorMap{ComplexSpace, 0, 2},V::AbstractTensorMap{ComplexSpace, 1, 2},nextψ::AbstractTensorMap{ComplexSpace, 1, 2})
    @tensor tempMPS[-1,-2,-3] ≔ Σ[1,-1]*nextψ[1,-2,-3]
    return [V,permute(tempMPS,(),(1,2,3))]
end

function LeftMerge(Σ::AbstractTensorMap{ComplexSpace, 0, 2},V::AbstractTensorMap{ComplexSpace, 1, 2},nextψ::AbstractTensorMap{ComplexSpace, 1, 2})
    @tensor tempMPS[-1,-2,-3] ≔ Σ[1,-3]*nextψ[1,-1,-2]
    return [permute(tempMPS,(),(1,2,3)),V]
end

function SVDMerge(
    Opri1::AbstractTensorMap{ComplexSpace, 1, 3},
    Opri2::AbstractTensorMap{ComplexSpace, 2, 2},
    svdOpri1::AbstractTensorMap{ComplexSpace, 2, 2},
    svdOpri2::AbstractTensorMap{ComplexSpace, 0, 2},
    )
    @tensor tempMPO[-1,-2,-3,-4] ≔ svdOpri2[1,-2]*Opri2[-1,1,-3,-4]

    return [svdOpri1,permute(tempMPO,(1,),(2,3,4))]
end

function SVDMerge(
    Opri1::AbstractTensorMap{ComplexSpace, 2, 2},
    Opri2::AbstractTensorMap{ComplexSpace, 1, 3},
    svdOpri1::AbstractTensorMap{ComplexSpace, 0, 2},
    svdOpri2::AbstractTensorMap{ComplexSpace, 2, 2},
    )
    @tensor tempMPO[-1,-2,-3,-4] ≔ Opri1[1,-1,-2,-3]*svdOpri1[1,-4]

    return [permute(tempMPO,(1,),(2,3,4)),svdOpri2]
end

function SVDMerge(
    Opri1::AbstractTensorMap{ComplexSpace, 2, 2},
    Opri2::AbstractTensorMap{ComplexSpace, 2, 2},
    svdOpri1::AbstractTensorMap{ComplexSpace, 0, 2},
    svdOpri2::AbstractTensorMap{ComplexSpace, 2, 2},
    )
    @tensor tempMPO[-1,-2,-3,-4] ≔ Opri1[1,-1,-2,-3]*svdOpri1[1,-4]

    return [permute(tempMPO,(1,),(2,3,4)),svdOpri2]
end

function SVDMerge(
    ψ1::AbstractTensorMap{ComplexSpace, 0, 3},
    ψ2::AbstractTensorMap{ComplexSpace, 1, 2},
    svdψ1::AbstractTensorMap{ComplexSpace, 1, 2},
    svdψ2::AbstractTensorMap{ComplexSpace, 0, 2},
    )
    @tensor tempMPS[-1,-2,-3] ≔ svdψ2[1,-1]*ψ2[1,-2,-3]

    return [svdψ1,permute(tempMPS,(),(1,2,3))]
end

function SVDMerge(
    ψ1::AbstractTensorMap{ComplexSpace, 1, 2},
    ψ2::AbstractTensorMap{ComplexSpace, 0, 3},
    svdψ1::AbstractTensorMap{ComplexSpace, 0, 2},
    svdψ2::AbstractTensorMap{ComplexSpace, 1, 2},
    )
    @tensor tempMPS[-1,-2,-3] ≔ svdψ1[1,-3]*ψ1[1,-1,-2]

    return [permute(tempMPS,(),(1,2,3)),svdψ2]
end

function EnvMerge(
    EnvL::AbstractTensorMap{ComplexSpace,1,1},
    Opr1::AbstractTensorMap{ComplexSpace,1,3},
    Opr2::AbstractTensorMap{ComplexSpace,2,2},
    EnvR::AbstractTensorMap{ComplexSpace,1,1}
    )

    @tensor tempMPO[-1,-2,-3,-4,-5,-6] ≔ EnvL[1,-3]*Opr1[-2,1,-4,2]*Opr2[-1,2,-5,3]*EnvR[3,-6]
    return permute(tempMPO,(1,2),(3,4,5,6))
end

function EnvMerge(
    EnvL::AbstractTensorMap{ComplexSpace,2,1},
    Opr11::AbstractTensorMap{ComplexSpace,1,3},
    Opr12::AbstractTensorMap{ComplexSpace,2,2},
    Opr21::AbstractTensorMap{ComplexSpace,1,3},
    Opr22::AbstractTensorMap{ComplexSpace,2,2},
    EnvR::AbstractTensorMap{ComplexSpace,2,1}
    )

    @tensor tempMPO[-1,-2,-3,-4,-5,-6] ≔ EnvL[1,2,-3]*Opr11[-2,1,3,4]*Opr21[3,2,-4,5]*Opr12[-1,4,6,7]*Opr22[6,5,-5,8]*EnvR[7,8,-6]
    return permute(tempMPO,(1,2),(3,4,5,6))
end

function EnvMerge(
    EnvL::AbstractTensorMap{ComplexSpace,1,1},
    Opr1::AbstractTensorMap{ComplexSpace,2,2},
    Opr2::AbstractTensorMap{ComplexSpace,1,3},
    EnvR::AbstractTensorMap{ComplexSpace,1,1}
    )

    @tensor tempMPO[-1,-2,-3,-4,-5,-6] ≔ EnvL[1,-3]*Opr1[2,-2,1,-4]*Opr2[-1,2,-5,3]*EnvR[3,-6]
    return permute(tempMPO,(1,2),(3,4,5,6))
end

function EnvMerge(
    EnvL::AbstractTensorMap{ComplexSpace,2,1},
    Opr11::AbstractTensorMap{ComplexSpace,2,2},
    Opr12::AbstractTensorMap{ComplexSpace,1,3},
    Opr21::AbstractTensorMap{ComplexSpace,2,2},
    Opr22::AbstractTensorMap{ComplexSpace,1,3},
    EnvR::AbstractTensorMap{ComplexSpace,2,1}
    )

    @tensor tempMPO[-1,-2,-3,-4,-5,-6] ≔ EnvL[7,8,-3]*Opr11[4,-2,7,6]*Opr21[5,6,8,-4]*Opr12[-1,4,3,1]*Opr22[3,5,-5,2]*EnvR[1,2,-6]
    return permute(tempMPO,(1,2),(3,4,5,6))
end

function EnvMerge(
    EnvL::AbstractTensorMap{ComplexSpace,1,1},
    ψ1::AbstractTensorMap{ComplexSpace,0,3},
    ψ2::AbstractTensorMap{ComplexSpace,1,2},
    EnvR::AbstractTensorMap{ComplexSpace,1,1}
    )

    @tensor tempMPS[-1,-2,-3,-4] ≔ EnvL[1,-1]*ψ1[1,-2,2]*ψ2[2,-3,3]*EnvR[3,-4]
    return permute(tempMPS,(),(1,2,3,4))
end

function EnvMerge(
    EnvL::AbstractTensorMap{ComplexSpace,1,1},
    ψ1::AbstractTensorMap{ComplexSpace,1,2},
    ψ2::AbstractTensorMap{ComplexSpace,0,3},
    EnvR::AbstractTensorMap{ComplexSpace,1,1}
    )

    @tensor tempMPS[-1,-2,-3,-4] ≔ EnvL[3,-1]*ψ1[2,3,-2]*ψ2[2,-3,1]*EnvR[1,-4]
    return permute(tempMPS,(),(1,2,3,4))
end


function EnvMerge(
    EnvL::AbstractTensorMap{ComplexSpace,2,1},
    localOpr::AbstractTensorMap{ComplexSpace,2,4},
    Opr1::AbstractTensorMap{ComplexSpace,1,3},
    Opr2::AbstractTensorMap{ComplexSpace,2,2},
    EnvR::AbstractTensorMap{ComplexSpace,2,1}
    )

    @tensor tempMPS[-1,-2,-3,-4,-5,-6] ≔ EnvL[1,2,-3]*localOpr[-1,-2,1,3,5,6]*Opr1[3,2,-4,4]*Opr2[5,4,-5,7]*EnvR[6,7,-6]
    return permute(tempMPS,(1,2),(3,4,5,6))
end

function EnvMerge(
    EnvL::AbstractTensorMap{ComplexSpace,2,1},
    localOpr::AbstractTensorMap{ComplexSpace,2,4},
    Opr1::AbstractTensorMap{ComplexSpace,2,2},
    Opr2::AbstractTensorMap{ComplexSpace,1,3},
    EnvR::AbstractTensorMap{ComplexSpace,2,1}
    )

    @tensor tempMPS[-1,-2,-3,-4,-5,-6] ≔ EnvL[6,7,-3]*Opr1[4,5,7,-4]*Opr2[3,4,-5,2]*localOpr[-1,-2,6,5,3,1]*EnvR[1,2,-6]
    return permute(tempMPS,(1,2),(3,4,5,6))
end

function EnvMerge(
    EnvL::AbstractTensorMap{ComplexSpace,2,1},
    localψ::AbstractTensorMap{ComplexSpace,0,4},
    Opr1::AbstractTensorMap{ComplexSpace,1,3},
    Opr2::AbstractTensorMap{ComplexSpace,2,2},
    EnvR::AbstractTensorMap{ComplexSpace,2,1}
    )
    @tensor tempMPS[-1,-2,-3,-4] ≔ EnvL[1,2,-1]*localψ[1,3,5,6]*Opr1[3,2,-2,4]*Opr2[5,4,-3,7]*EnvR[6,7,-4]
    return permute(tempMPS,(),(1,2,3,4))
end

function EnvMerge(
    EnvL::AbstractTensorMap{ComplexSpace,2,1},
    localψ::AbstractTensorMap{ComplexSpace,0,4},
    Opr1::AbstractTensorMap{ComplexSpace,2,2},
    Opr2::AbstractTensorMap{ComplexSpace,1,3},
    EnvR::AbstractTensorMap{ComplexSpace,2,1}
    )

    @tensor tempMPS[-1,-2,-3,-4,] ≔ EnvL[6,7,-1]*Opr1[4,5,7,-2]*Opr2[3,4,-3,2]*localψ[6,5,3,1]*EnvR[1,2,-4]
    return permute(tempMPS,(),(1,2,3,4))
end

function EnvMerge(
    EnvL::AbstractTensorMap{ComplexSpace,2,1},
    localOpr::AbstractTensorMap{ComplexSpace,1,3},
    Opr::AbstractTensorMap{ComplexSpace,1,3},
    EnvR::AbstractTensorMap{ComplexSpace,2,1}
    )

    @tensor tempMPS[-1,-2,-3,-4] ≔ EnvL[1,2,-2]*localOpr[-1,1,3,4]*Opr[3,2,-3,5]*EnvR[4,5,-4]
    return permute(tempMPS,(1,),(2,3,4))
end

########### NON - normalized ###########

function EnvMerge(
    EnvL::AbstractTensorMap{ComplexSpace,2,1},
    localψ::AbstractTensorMap{ComplexSpace,0,4},
    Opr1::AbstractTensorMap{ComplexSpace,2,2},
    Opr2::AbstractTensorMap{ComplexSpace,2,2},
    EnvR::AbstractTensorMap{ComplexSpace,1,2}
    )

    @tensor tempMPS[-1,-2,-3,-4] ≔ EnvL[1,2,-1]*localψ[1,3,5,6]*Opr1[4,3,2,-2]*Opr2[7,5,4,-3]*EnvR[6,7,-4]
    return permute(tempMPS,(),(1,2,3,4))
end

function EnvMerge(
    EnvL::AbstractTensorMap{ComplexSpace,2,1},
    localψ::AbstractTensorMap{ComplexSpace,0,3},
    Opr::AbstractTensorMap{ComplexSpace,2,2},
    EnvR::AbstractTensorMap{ComplexSpace,1,2}
    )

    @tensor tempMPS[-1,-2,-3] ≔ EnvL[1,2,-1]*localψ[1,3,4]*Opr[5,3,2,-2]*EnvR[4,5,-3]
    return permute(tempMPS,(),(1,2,3))
end

