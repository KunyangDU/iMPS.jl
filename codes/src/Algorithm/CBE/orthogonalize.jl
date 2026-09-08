
# 逐元素归并 per-worker 私有输出向量（键索引 → 累加张量），供 threaded_reduce! 的 combine! 使用
function _merge_envs!(x::Vector{Any}, y::Vector{Any})
    for i in eachindex(x)
        if x[i] === nothing
            x[i] = y[i]
        elseif y[i] !== nothing
            axpy!(1, y[i], x[i])
        end
    end
    return x
end

function orthogonalize!(B::T,H::SparseMPOTensor,B′::T′,EnvR::SparseRightEnvironmentTensor) where T <: Union{DenseMPOTensor{4}, MPSTensor{3}} where T′ <: Union{AdjointMPOTensor{4}, AdjointMPSTensor{3}}
    # 算子优先（元任务 = _validind 的每个算符 j）：预求和右环境 → 一次缩并 → 按左键散射。
    # 旧 bond-first 版本同一算符的每个左键都重复缩并一次；这里每算符只缩并一次，任务数 = 算符数（≫ 键数）。
    # 散射用 per-worker 私有输出向量 + 末步归并，无嵌套循环、无 barrier、不物化整份 tmpC。
    validind = _validind(H)
    accs = [Vector{Any}(nothing, length(H.left.fwd)) for _ in 1:get_nworker()]
    merged = threaded_reduce!(eachindex(validind), accs; combine! = _merge_envs!) do k, acc, _
        l_inds, j, r_inds, wl, wr = validind[k]
        weighted_env = _wsum(EnvR, r_inds, wr)
        C = rorth!(contract(H[j], B′, weighted_env), B)
        for (idx, i) in enumerate(l_inds)
            acc[i] = axpy!(wl[idx], C, acc[i])
        end
        acc
    end
    return SparseRightEnvironmentTensor(convert(Vector{RightCompositeEnvironmentTensor}, merged))
end

function orthogonalize!(A::T,H::SparseMPOTensor,A′::T′,EnvL::SparseLeftEnvironmentTensor) where T <: Union{DenseMPOTensor{4},MPSTensor{3}} where T′ <: Union{AdjointMPOTensor{4}, AdjointMPSTensor{3}}
    # 算子优先（pushright 镜像）：预求和左环境 → 一次缩并 → 按右键散射，同样每算符只缩并一次。
    validind = _validind(H)
    accs = [Vector{Any}(nothing, length(H.right.rev)) for _ in 1:get_nworker()]
    merged = threaded_reduce!(eachindex(validind), accs; combine! = _merge_envs!) do k, acc, _
        l_inds, j, r_inds, wl, wr = validind[k]
        weighted_env = _wsum(EnvL, l_inds, wl)
        C = lorth!(contract(weighted_env, H[j], A′), A)
        for (idx, i) in enumerate(r_inds)
            acc[i] = axpy!(wr[idx], C, acc[i])
        end
        acc
    end
    return SparseLeftEnvironmentTensor(convert(Vector{LeftCompositeEnvironmentTensor}, merged))
end

# function orthogonalize!(A::Union{DenseMPOTensor{4},MPSTensor{3}},A′::Union{DenseMPOTensor{4},MPSTensor{3}},Env::Union{DenseLeftEnvironmentTensor,DenseRightEnvironmentTensor})
#     tmp = contract(Env.A,A)
#     Envorth = tmp - contract(tmp,A′)
#     return Envorth
# end

# function orthogonalize!(Q::T,A::T,direction::AbstractDirection;tol::Number=1e-16) where T <: Union{MPSTensor{3},DenseMPOTensor{4},AdjointMPOTensor{4}}
#     norm(Q) ≈ 0 && return Q
#     ϵ = norm(_cbeinner(Q,A,direction)) / norm(Q)
#     for _ in 1:10
#         ϵ = _cbeorth!(Q,A,direction) / norm(Q)
#         ϵ < tol && break
#     end
#     @assert ϵ < tol ϵ
#     return Q
# end

orthogonalize!(A::T₁,H::T,A′::T₂,EnvL::DenseLeftEnvironmentTensor) where {T <: Union{DenseMPOTensor{4},AdjointMPOTensor{4}}, T₁ <: Union{DenseMPOTensor{4}, MPSTensor{3}}, T₂ <: Union{AdjointMPOTensor{4},AdjointMPSTensor{3}}} = _orth_sub!(contract(EnvL.A, H, A′), A)
orthogonalize!(B::T₁,H::T,B′::T₂,EnvR::DenseRightEnvironmentTensor) where {T <: Union{DenseMPOTensor{4},AdjointMPOTensor{4}}, T₁ <: Union{DenseMPOTensor{4}, MPSTensor{3}}, T₂ <: Union{AdjointMPOTensor{4},AdjointMPSTensor{3}}} = _orth_sub!(contract(H, B′, EnvR.A), B)

