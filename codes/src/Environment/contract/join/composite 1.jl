function contract(El::LeftCompositeEnvironmentTensor{1, 3, 3, 1}, Er::RightEnvironmentTensor{2})
    @tensor tmp[-1;-2 -3] ≔ El.A[1,-3,-2] * Er.A[-1,1]
    return AdjointMPSTensor(tmp)
end

function contract(El::LeftEnvironmentTensor{2}, Er::RightCompositeEnvironmentTensor{2, 3, 3, 1})
    @tensor tmp[-1;-2 -3] ≔ El.A[1,-2] *Er.A[-1,1,-3]
    return AdjointMPSTensor(tmp)
end

function contract(El::LeftCompositeEnvironmentTensor{1, 4, 3, 1}, Er::RightEnvironmentTensor{3})
    @tensor tmp[-1;-2 -3] ≔ El.A[1,2,-3,-2] * Er.A[-1,2,1]
    return AdjointMPSTensor(tmp)
end

function contract(El::LeftEnvironmentTensor{3}, Er::RightCompositeEnvironmentTensor{2, 4, 3, 1})
    @tensor tmp[-1;-2 -3] ≔ El.A[1,2,-2] * Er.A[-1,2,1,-3]
    return AdjointMPSTensor(tmp)
end

contract(El::LeftCompositeEnvironmentTensor{2, 4, 3, 1}, Er::RightCompositeEnvironmentTensor{2, 4, 3, 1}) = AdjointCompositeMPOTensor(@tensor tmp[-1 -2 -3;-4 -5 -6] ≔ El.A[1,-3,-5,-6] * Er.A[-1,-2,1,-4])
