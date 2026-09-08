function lorth!(x::LeftCompositeEnvironmentTensor{1, 3, 3, 1}, A::MPSTensor{3})
    @tensor x.A[-1;-4 -5] -= x.A[-1,2,1] * A.A[1,2,3] * A.A'[3,-5,-4]
    return x
end

function lorth!(x::LeftCompositeEnvironmentTensor{1, 4, 3, 1}, A::MPSTensor{3})
    @tensor x.A[-1;-3 -4 -5] -= x.A[-1,-3,2,1] * A.A[1,2,3] * A.A'[3,-5,-4]
    return x
end

function lorth!(x::LeftCompositeEnvironmentTensor{2, 4, 3, 1}, A::DenseMPOTensor{4})
    @tensor x.A[-1 -2;-4 -5] -= x.A[-1,3,2,1] * A.A[2,1,4,3] * A.A'[4,-2,-4,-5]
    return x
end

function rorth!(x::RightCompositeEnvironmentTensor{2, 3, 3, 1}, A::MPSTensor{3})
    @tensor x.A[-1;-4 -5] -= x.A[1,-4,2] * A.A[3,2,1] * A.A'[-1,3,-5]
    return x
end

function rorth!(x::RightCompositeEnvironmentTensor{2, 4, 3, 1}, A::MPSTensor{3})
    @tensor x.A[-1 -2;-4 -5] -= x.A[1,-2,-4,2] * A.A[3,2,1] * A.A'[-1,3,-5]
    return x
end

function rorth!(x::RightCompositeEnvironmentTensor{2, 4, 3, 1}, A::DenseMPOTensor{4})
    @tensor x.A[-1 -2;-4 -5] -= x.A[1,3,-4,2] * A.A[2,4,1,3] * A.A'[-1,-2,-5,4]
    return x
end