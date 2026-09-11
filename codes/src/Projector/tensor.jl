function lproj!(obj::CompositeMPSTensor{2,4}, A::MPSTensor{3})
    @tensor obj.A[-1 -2 -3;-4] = obj.A[1,2,-3,-4] * A.A'[3,1,2] * A.A[-1,-2,3]
    return obj
end

function lorth!(obj::CompositeMPSTensor{2,4}, A::MPSTensor{3})
    @tensor obj.A[-1 -2 -3;-4] -= obj.A[1,2,-3,-4] * A.A'[3,1,2] * A.A[-1,-2,3]
    return obj
end

function rproj!(obj::CompositeMPSTensor{2,4}, A::MPSTensor{3})
    @tensor obj.A[-1 -2 -3;-4] = obj.A[-1,-2,2,1] * A.A'[1,3,2] * A.A[3,-3,-4]
    return obj
end

function rorth!(obj::CompositeMPSTensor{2,4}, A::MPSTensor{3})
    @tensor obj.A[-1 -2 -3;-4] -= obj.A[-1,-2,2,1] * A.A'[1,3,2] * A.A[3,-3,-4]
    return obj
end

function lproj!(obj::CompositeMPOTensor{2,6}, A::DenseMPOTensor{4})
    @tensor obj.A[-1 -2 -3;-4 -5 -6] = obj.A[-1,2,1,-4,-5,3] * A.A'[4,3,2,1] * A.A[-2,-3,4,-6]
    return obj
end

function rproj!(obj::CompositeMPOTensor{2,6}, A::DenseMPOTensor{4})
    @tensor obj.A[-1 -2 -3;-4 -5 -6] = obj.A[2,-2,-3,1,3,-6] * A.A'[1,3,2,4] * A.A[-1,4,-4,-5]
    return obj
end

function lorth!(obj::CompositeMPOTensor{2,6}, A::DenseMPOTensor{4})
    @tensor obj.A[-1 -2 -3;-4 -5 -6] -= obj.A[-1,2,1,-4,-5,3] * A.A'[4,3,2,1] * A.A[-2,-3,4,-6]
    return obj
end

function rorth!(obj::CompositeMPOTensor{2,6}, A::DenseMPOTensor{4})
    @tensor obj.A[-1 -2 -3;-4 -5 -6] -= obj.A[2,-2,-3,1,3,-6] * A.A'[1,3,2,4] * A.A[-1,4,-4,-5]
    return obj
end

function lproj(obj::CompositeMPOTensor{2,6}, A::DenseMPOTensor{4})
    @tensor tmp[-1 -2 -3;-4 -5 -6] ≔ obj.A[-1,2,1,-4,-5,3] * A.A'[4,3,2,1] * A.A[-2,-3,4,-6]
    return CompositeMPOTensor(tmp)
end

function rproj(obj::CompositeMPOTensor{2,6}, A::DenseMPOTensor{4})
    @tensor tmp[-1 -2 -3;-4 -5 -6] ≔ obj.A[2,-2,-3,1,3,-6] * A.A'[1,3,2,4] * A.A[-1,4,-4,-5]
    return CompositeMPOTensor(tmp)
end

function lorth(obj::CompositeMPOTensor{2,6}, A::DenseMPOTensor{4})
    @tensor tmp[-1 -2 -3;-4 -5 -6] ≔ obj.A[-1,2,1,-4,-5,3] * A.A'[4,3,2,1] * A.A[-2,-3,4,-6]
    return CompositeMPOTensor(obj.A - tmp)
end

function rorth(obj::CompositeMPOTensor{2,6}, A::DenseMPOTensor{4})
    @tensor tmp[-1 -2 -3;-4 -5 -6] ≔ obj.A[2,-2,-3,1,3,-6] * A.A'[1,3,2,4] * A.A[-1,4,-4,-5]
    return CompositeMPOTensor(obj.A - tmp)
end

