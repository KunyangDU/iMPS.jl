function contract(El::LeftEnvironmentTensor{3}, h::DenseMPOTensor{4}, A::AdjointMPSTensor{3})
    @tensor tmp[-1;-3 -4 -5] ≔ El.A[1,3,-5] * h.A[2,3,-3,-4] * A.A[-1,1,2]
    return LeftCompositeEnvironmentTensor(tmp,3,1)
end

function contract(h::DenseMPOTensor{4}, A::AdjointMPSTensor{3}, Er::RightEnvironmentTensor{3})
    @tensor tmp[-1 -2;-4 -5] ≔ h.A[3,-2,2,-5] * A.A[1,-4,3] * Er.A[-1,2,1]
    return RightCompositeEnvironmentTensor(tmp,3,1)
end