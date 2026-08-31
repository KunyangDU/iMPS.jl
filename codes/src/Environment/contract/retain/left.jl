contract(EnvL::LeftCompositeEnvironmentTensor{2, 3, 3, 3}, A::AdjointMPSTensor{3}) = LeftEnvironmentTensor(@tensor tmp[-1;-2] ≔ EnvL.A[1,2,-2] * A.A[-1,1,2] )
contract(EnvL::LeftCompositeEnvironmentTensor{2, 4, 3, 3}, A::AdjointMPSTensor{3}) = LeftEnvironmentTensor(@tensor tmp[-1;-2 -3] ≔ EnvL.A[1,2,-2,-3] * A.A[-1,1,2] )

contract(EnvL::LeftCompositeEnvironmentTensor{2, 4, 3, 3}, A::AdjointMPOTensor{4}) = LeftEnvironmentTensor(@tensor tmp[-1;-2] ≔ EnvL.A[1,2,-2,3] * A.A[-1,3,2,1] )
contract(EnvL::LeftCompositeEnvironmentTensor{2, 5, 3, 3}, A::AdjointMPOTensor{4}) = LeftEnvironmentTensor(@tensor tmp[-1;-2 -3] ≔ EnvL.A[1,2,-2,-3,3] * A.A[-1,3,2,1] )

function contract(El::LeftCompositeEnvironmentTensor{1, 3, 3, 1}, A::MPSTensor{3})
    @tensor tmp[-1;-2] ≔ El.A[-1,2,1] * A.A[1,2,-2]
    return LeftEnvironmentTensor(tmp)
end

function contract(El::LeftCompositeEnvironmentTensor{1, 4, 3, 1}, A::MPSTensor{3})
    @tensor tmp[-1;-2 -3] ≔ El.A[-1,-2,2,1] * A.A[1,2,-3]
    return LeftEnvironmentTensor(tmp)
end
