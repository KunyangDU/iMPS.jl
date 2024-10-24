
#################### PUSH ####################

function PushLeft(EnvR::AbstractTensorMap{ComplexSpace,1,2},ψi::AbstractTensorMap{ComplexSpace,1,2},Hi::AbstractTensorMap{ComplexSpace,2,2})
    @tensor EnvRR[-1,-2,-3] ≔ ψi[-1,4,1]*Hi[3,4,-2,5]*ψi'[5,2,-3]*EnvR[1,3,2]
    return permute(EnvRR,(1,),(2,3))
end

function PushLeft(EnvR::AbstractTensorMap{ComplexSpace,1,2},ψi1::AbstractTensorMap{ComplexSpace,1,2},Hi::AbstractTensorMap{ComplexSpace,2,2},ψi2::AbstractTensorMap{ComplexSpace,1,2})
    @tensor EnvRR[-1,-2,-3] ≔ ψi1[-1,4,1]*Hi[3,4,-2,5]*ψi2'[5,2,-3]*EnvR[1,3,2]
    return permute(EnvRR,(1,),(2,3))
end

#################### normalized ####################

function PushLeft(EnvR::AbstractTensorMap{ComplexSpace,1,1},ψi1::AbstractTensorMap{ComplexSpace,1,2},ψi2::AbstractTensorMap{ComplexSpace,1,2})
    @tensor EnvRR[-1,-2] ≔ ψi1[-1,3,1] * ψi2'[3,2,-2] * EnvR[1,2]
    return permute(EnvRR,(1,),(2,))
end

function PushLeft(EnvR::AbstractTensorMap{ComplexSpace,2,1},ψi::AbstractTensorMap{ComplexSpace,1,2},Opri::AbstractTensorMap{ComplexSpace,2,2})
    @tensor EnvRR[-1,-2,-3] ≔ ψi[-1,3,1] * Opri[3,-2,5,2] * ψi'[5,4,-3] * EnvR[1,2,4]
    return permute(EnvRR,(1,2),(3,))
end

function PushLeft(EnvR::AbstractTensorMap{ComplexSpace,2,1},ψi1::AbstractTensorMap{ComplexSpace,1,2},Opri::AbstractTensorMap{ComplexSpace,2,2},ψi2::AbstractTensorMap{ComplexSpace,1,2})
    @tensor EnvRR[-1,-2,-3] ≔ ψi1[-1,3,1] * Opri[3,-2,5,2] * ψi2'[5,4,-3] * EnvR[1,2,4]
    return permute(EnvRR,(1,2),(3,))
end

function PushLeft(EnvR::AbstractTensorMap{ComplexSpace,1,1},Opri1::AbstractTensorMap{ComplexSpace,2,2},Opri2::AbstractTensorMap{ComplexSpace,2,2})
    @tensor EnvRR[-1,-2] ≔ Opri1[4,-1,3,1] * Opri2'[3,2,4,-2] * EnvR[1,2]
    return permute(EnvRR,(1,),(2,))
end

function PushLeft(EnvR::AbstractTensorMap{ComplexSpace,2,1},Opri1::AbstractTensorMap{ComplexSpace,2,2},Opri2::AbstractTensorMap{ComplexSpace,2,2},Opri3::AbstractTensorMap{ComplexSpace,2,2})
    @tensor EnvRR[-1,-2,-3] ≔ Opri1[5,-1,3,1] * Opri2[3,-2,6,2] * Opri3'[6,4,5,-3] * EnvR[1,2,4]
    return permute(EnvRR,(1,2),(3,))
end

function PushLeft(EnvR::AbstractTensorMap{ComplexSpace,1,0},Opri::AbstractTensorMap{ComplexSpace,2,2})
    @tensor EnvRR[-1] ≔ Opri[2,-1,2,1] * EnvR[1]
    return permute(EnvRR,(1,),())
end