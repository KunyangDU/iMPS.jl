function contract(El::LeftEnvironmentTensor{2}, h::LocalOperator{1,1}, A::AdjointMPSTensor{3})
    @tensor tmp[-1;-4 -5] ≔ El.A[1,-5] * h.A[2,-4] * A.A[-1,1,2]
    return LeftCompositeEnvironmentTensor(tmp,3,1)
end

function contract(El::LeftEnvironmentTensor{2}, ::IdentityOperator{1}, A::AdjointMPSTensor{3})
    @tensor tmp[-1;-4 -5] ≔ El.A[1,-5] * A.A[-1,1,-4]
    return LeftCompositeEnvironmentTensor(tmp,3,1)
end

function contract(h::LocalOperator{1, 1}, A::AdjointMPSTensor{3}, Er::RightEnvironmentTensor{2})
    @tensor tmp[-1;-4 -5] ≔ h.A[2,-5] * A.A[1,-4,2] * Er.A[-1,1]
    return RightCompositeEnvironmentTensor(tmp,3,1)
end

function contract(::IdentityOperator{1}, A::AdjointMPSTensor{3}, Er::RightEnvironmentTensor{2})
    @tensor tmp[-1;-4 -5] ≔ A.A[1,-4,-5] * Er.A[-1,1]
    return RightCompositeEnvironmentTensor(tmp,3,1)
end
