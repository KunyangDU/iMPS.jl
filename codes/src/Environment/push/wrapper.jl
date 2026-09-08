function pushleft!(env::Environment{R}) where R
    @assert 1 ≤ env.center[1] ≤ env.center[2] ≤ env.L

    env.envs[env.center[2]] = pushleft(env.layer..., env.envs[env.center[2] + 1], env.center[2])

    env.center[2] -= 1
    ( env.center[1] > env.center[2] ) && ( env.center[1] -= 1 )
end

function pushright!(env::Environment{R}) where R
    @assert 1 ≤ env.center[1] ≤ env.center[2] ≤ env.L

    env.envs[env.center[1] + 1] = pushright(env.layer..., env.envs[env.center[1]], env.center[1])

    env.center[1] += 1
    ( env.center[1] > env.center[2] ) && ( env.center[2] += 1 )
end

pushleft(A::DenseMPS, mpo::SparseMPO, B::T, EnvR::SparseRightEnvironmentTensor{1}, i::Int64) where T <: Union{AdjointMPS,RefMPS} = pushleft(A[i],mpo[i],B[i],EnvR)
pushright(A::DenseMPS, mpo::SparseMPO, B::T, EnvL::SparseLeftEnvironmentTensor{1}, i::Int64) where T <: Union{AdjointMPS,RefMPS} = pushright(A[i],mpo[i],B[i],EnvL)

pushleft(A::DenseMPO, B::SparseMPO, C::T, EnvR::SparseRightEnvironmentTensor, i::Int64) where T <: Union{AdjointMPO,RefMPO} = pushleft(A[i],B[i],C[i],EnvR)
pushright(A::DenseMPO, B::SparseMPO, C::T, EnvL::SparseLeftEnvironmentTensor, i::Int64) where T <: Union{AdjointMPO,RefMPO} = pushright(A[i],B[i],C[i],EnvL)

pushleft(A::DenseMPO, B::DenseMPO, C::T, EnvR::DenseRightEnvironmentTensor{3}, site::Int64) where T <: Union{AdjointMPO,RefMPO} = DenseRightEnvironmentTensor(contract(A[site], B[site], C[site], EnvR.A))
pushright(A::DenseMPO, B::DenseMPO, C::T, EnvL::DenseLeftEnvironmentTensor{3}, site::Int64) where T <: Union{AdjointMPO,RefMPO} = DenseLeftEnvironmentTensor(contract(A[site], B[site], C[site], EnvL.A))

pushleft(A::DenseMPS, B::DenseMPO, C::T₃, EnvR::DenseRightEnvironmentTensor{3}, site::Int64) where T₃ <: Union{AdjointMPS,RefMPS} = DenseRightEnvironmentTensor(contract(A[site], B[site], C[site], EnvR.A))
pushright(A::DenseMPS, B::DenseMPO, C::T₃, EnvL::DenseLeftEnvironmentTensor{3}, site::Int64) where T₃ <: Union{AdjointMPS,RefMPS} = DenseLeftEnvironmentTensor(contract(A[site], B[site], C[site], EnvL.A))

pushleft(A::DenseMPO, B::Union{AdjointMPO,RefMPO}, C::AdjointMPO, EnvR::DenseRightEnvironmentTensor{3}, site::Int64) = DenseRightEnvironmentTensor(contract(A[site], B[site], C[site], EnvR.A))
pushright(A::DenseMPO, B::Union{AdjointMPO,RefMPO}, C::AdjointMPO, EnvL::DenseLeftEnvironmentTensor{3}, site::Int64) = DenseLeftEnvironmentTensor(contract(A[site], B[site], C[site], EnvL.A))

# layer 2 - dense

pushleft(A::DenseMPS, B::Union{AdjointMPS,RefMPS}, EnvR::DenseRightEnvironmentTensor{2}, site::Int64) = DenseRightEnvironmentTensor(contract(A[site], B[site], EnvR.A))
pushright(A::DenseMPS, B::Union{AdjointMPS,RefMPS}, EnvL::DenseLeftEnvironmentTensor{2}, site::Int64) = DenseLeftEnvironmentTensor(contract(A[site], B[site], EnvL.A))

pushleft(A::DenseMPO, B::Union{AdjointMPO,RefMPO}, EnvR::DenseRightEnvironmentTensor{2}, site::Int64) = DenseRightEnvironmentTensor(contract(A[site], B[site], EnvR.A))
pushright(A::DenseMPO, B::Union{AdjointMPO,RefMPO}, EnvL::DenseLeftEnvironmentTensor{2}, site::Int64) = DenseLeftEnvironmentTensor(contract(A[site], B[site], EnvL.A))


