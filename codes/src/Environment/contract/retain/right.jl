contract(EnvR::RightCompositeEnvironmentTensor{1, 3, 3, 3}, A::AdjointMPSTensor{3}) = RightEnvironmentTensor(@tensor tmp[-1;-2] ≔ EnvR.A[-1,2,1] * A.A[1,-2,2] )
contract(EnvR::RightCompositeEnvironmentTensor{1, 4, 3, 3}, A::AdjointMPSTensor{3}) = RightEnvironmentTensor(@tensor tmp[-1,-2;-3] ≔ EnvR.A[-1,-2,2,1] * A.A[1,-3,2])

contract(EnvR::RightCompositeEnvironmentTensor{2, 4, 3, 3}, A::AdjointMPOTensor{4}) = RightEnvironmentTensor(@tensor tmp[-1;-2] ≔ EnvR.A[-1,2,1,3] * A.A[1,3,2,-2] )
contract(EnvR::RightCompositeEnvironmentTensor{2, 5, 3, 3}, A::AdjointMPOTensor{4}) = RightEnvironmentTensor(@tensor tmp[-1 -2;-3] ≔ EnvR.A[-1,-2,2,1,3] * A.A[1,3,2,-3] )

function contract(Er::RightCompositeEnvironmentTensor{2, 3, 3, 1}, A::MPSTensor{3})
    @tensor tmp[-1;-2] ≔ A.A[-1,2,1] * Er.A[1,-2,2]
    return RightEnvironmentTensor(tmp)
end

function contract(Er::RightCompositeEnvironmentTensor{2, 4, 3, 1}, A::MPSTensor{3})
    @tensor tmp[-1 -2;-3] ≔ A.A[-1,2,1] * Er.A[1,-2,-3,2]
    return RightEnvironmentTensor(tmp)
end
