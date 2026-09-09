# ====================== wrapper 通用函数补全 ======================
# 覆盖三类对象:
#   1. AbstractTensorWrapper (MPSTensor / DenseMPOTensor / Composite 张量 / Environment 张量 ...)
#      —— 按 VectorInterface + LinearAlgebra 协议补全, 直接扩展 TensorKit (= VI) 的泛函;
#      限定写法 (TensorKit.xxx / LinearAlgebra.xxx) 保证 KrylovKit 内部派发命中, 不在 Main 遮蔽。
#   2. MP 层容器 (DenseMPS / AdjointMPS / DenseMPO / AdjointMPO / RefMPS / RefMPO / SparseMPO)
#      —— 索引 / 长度 / 归一化 (normalize! 归一化同时返回 norm, 区别于 LinearAlgebra 约定)。
#   3. Sparse 环境张量与 CompositeMPOTensor 的专用就地代数。
#
# 注意: 本文件在 using TensorKit 之后 include (见 TenetKit.jl), 未限定定义会遮蔽
# TensorKit re-export 的泛函 —— 所有跨包泛函一律带命名空间前缀。

# ---------- 标量与空间 ----------

TensorKit.scalartype(A::AbstractTensorWrapper) = TensorKit.scalartype(A.A)
TensorKit.scalartype(::Type{W}) where {W<:AbstractTensorWrapper} = scalartype(fieldtype(W, :A))

Base.eltype(A::AbstractTensorWrapper) = eltype(A.A)
TensorKit.space(A::AbstractTensorWrapper) = space(A.A)
TensorKit.space(A::AbstractLocalOperator) = space(A.A)
TensorKit.space(A::AbstractTensorWrapper, i::Int64) = space(A.A, i)
TensorKit.space(A::AbstractLocalOperator, i::Int64) = space(A.A, i)
TensorKit.codomain(A::AbstractTensorWrapper) = codomain(A.A)
TensorKit.domain(A::AbstractTensorWrapper) = domain(A.A)
TensorKit.dims(A::AbstractTensorWrapper) = dims(A.A)
TensorKit.numind(A::AbstractTensorWrapper) = numind(A.A)

# ---------- 零向量 / 转换 / 拷贝 ----------

# zerovector 只服务自同态场景 (输出空间 == 输入空间); 非自同态用 actionb。
function TensorKit.zerovector(A::T, ::Type{S}) where {S<:Number, T<:AbstractTensorWrapper}
    return convert(T, TensorKit.zerovector(A.A, S))
end
function TensorKit.zerovector!(A::AbstractTensorWrapper)
    TensorKit.zerovector!(A.A)
    return A
end
TensorKit.zerovector!!(A::AbstractTensorWrapper) = TensorKit.zerovector!(A)

Base.similar(A::AbstractTensorWrapper, ::Type{S}) where {S<:Number} = zerovector(A, S)
Base.convert(::Type{T}, A::AbstractTensorMap) where {T<:AbstractTensorWrapper} = T(A)
Base.copy(A::T) where {T<:AbstractTensorWrapper} = T(copy(A.A))
Base.iterate(t::AbstractTensorWrapper) = (t.A, nothing)
Base.iterate(::AbstractTensorWrapper, ::Nothing) = nothing

# ---------- 内积 / 范数 / 归一化 ----------

TensorKit.inner(A::T, B::T) where {T<:AbstractTensorWrapper} = inner(A.A, B.A)
TensorKit.norm(A::AbstractTensorWrapper) = norm(A.A)

# 本库约定: normalize! 就地归一化并返回 norm (与 LinearAlgebra 的返回对象不同)。
function TensorKit.LinearAlgebra.normalize!(obj::AbstractTensorWrapper)
    n = norm(obj.A)
    obj.A = obj.A / n
    return n
end
TensorKit.LinearAlgebra.normalize(A::T) where {T<:AbstractTensorWrapper} = T(normalize(A.A))

# ---------- 缩放 (scale 族) ----------

TensorKit.scale(A::T, α::Number) where {T<:AbstractTensorWrapper} = T(α * A.A)
TensorKit.scale!(A::AbstractTensorWrapper, α::Number) = TensorKit.LinearAlgebra.rmul!(A, α)

# 提升路径: α 与现有标量类型不兼容时分配新对象, 否则就地。
function TensorKit.scale!!(A::AbstractTensorWrapper, α::Number)
    S = promote_type(scalartype(A.A), typeof(α))
    return S <: scalartype(A.A) ? scale!(A, α) : scale(A, α)
end

function TensorKit.LinearAlgebra.rmul!(A::AbstractTensorWrapper, α::Number)
    TensorKit.LinearAlgebra.rmul!(A.A, α)
    return A
end
function TensorKit.LinearAlgebra.mul!(A::T, B::T, α::Number) where {T<:AbstractTensorWrapper}
    TensorKit.LinearAlgebra.mul!(A.A, B.A, α)
    return A
end

# ---------- 累加 (axpy / axpby / add 族) ----------
# axpby!(α, A, β, B) = B ← α*A + β*B (与 TensorKit/LinearAlgebra.axpby! 同序)。
# 就地路径要求 B 的标量类型够宽; 需要类型提升的场合用 add!!。

function TensorKit.LinearAlgebra.axpy!(α::Number, A::T, B::T) where {T<:AbstractTensorWrapper}
    TensorKit.axpy!(α, A.A, B.A)   # 零分配：要求 B.A 类型够宽（需提升用 add!!）
    return B
end

function TensorKit.LinearAlgebra.axpby!(α::Number, A::AbstractTensorWrapper, β::Number, B::AbstractTensorWrapper)
    TensorKit.axpby!(α, A.A, β, B.A)   # 零分配：要求 B.A 类型够宽（需提升用 add!!）
    return B
end
TensorKit.LinearAlgebra.axpby!(α::Number, A::AbstractTensorWrapper, ::Number, ::Nothing) = α * A
TensorKit.LinearAlgebra.axpby!(::Number, ::Nothing, β::Number, A::AbstractTensorWrapper) = TensorKit.LinearAlgebra.rmul!(A, β)
TensorKit.LinearAlgebra.axpy!(α::Number, A::AbstractTensorWrapper, ::Nothing) = α * A
TensorKit.LinearAlgebra.axpy!(::Number, ::Nothing, B::AbstractTensorWrapper) = B

TensorKit.add!(A::AbstractTensorWrapper, B::AbstractTensorWrapper) = TensorKit.LinearAlgebra.axpy!(true, B, A)
TensorKit.add!(A::AbstractTensorWrapper, ::Nothing) = A
TensorKit.add!(::Nothing, A::AbstractTensorWrapper) = A

# add!!(A, B, β, α) = A ← β*B + α*A; 语义与 VectorInterface.add!!(y, x, α, β) = α*x + β*y 一致。
# 一律外积(out-of-place)返回新对象: 就地路径会把 action(cache) 返回的持久累加器 c.acc[1] 与
# KrylovKit 的 Lanczos 基向量别名, 下一轮 action 的 zerovector! 逐次覆写基向量 → DMRG 能量错误。
function TensorKit.add!!(A::AbstractTensorWrapper,
                         B::AbstractTensorWrapper,
                         β::Number = one(scalartype(B)),
                         α::Number = one(scalartype(A)))
    return α * A + β * B
end
TensorKit.add!!(A::AbstractTensorWrapper, ::Nothing, α::Number = one(scalartype(A))) = α * A
TensorKit.add!!(::Nothing, B::AbstractTensorWrapper, β::Number = one(scalartype(B))) = β * B

# ---------- 初等运算 ----------

Base.:+(A::T, B::T) where {T<:AbstractTensorWrapper} = T(A.A + B.A)
Base.:+(::Nothing, B::AbstractTensorWrapper) = B
Base.:+(A::AbstractTensorWrapper, ::Nothing) = A
Base.:-(A::T, B::T) where {T<:AbstractTensorWrapper} = T(A.A - B.A)
Base.:-(A::T) where {T<:AbstractTensorWrapper} = T(-A.A)
Base.:*(A::T, B::T) where {T<:AbstractTensorWrapper} = T(A.A * B.A)
Base.:*(A::Number, B::T) where {T<:AbstractTensorWrapper} = T(A * B.A)
Base.:*(B::T, A::Number) where {T<:AbstractTensorWrapper} = T(A * B.A)
Base.:/(A::T, B::Number) where {T<:AbstractTensorWrapper} = (1 / B) * A

Base.isapprox(A::AbstractTensorWrapper, B::AbstractTensorWrapper) = isapprox(A.A, B.A)

# randn/rand 族: 直接复用 TensorKit 的 (Type, codomain, domain) 构造器。
Base.randn(A::T) where {T<:AbstractTensorWrapper} = T(randn(eltype(A), codomain(A), domain(A)))
Base.rand(A::T) where {T<:AbstractTensorWrapper} = T(rand(eltype(A), codomain(A), domain(A)))

# ---------- MP 层容器: 谓词 / 索引 / 长度 ----------

issparse(::T) where {T<:Union{DenseMPS,AdjointMPS,DenseMPO,AdjointMPO}} = false
issparse(::SparseMPO) = true
issparse(::SparseMPOTensor) = true

_isdisk(obj::T) where {T<:Union{DenseMPS,AdjointMPS,DenseMPO,AdjointMPO}} = obj.isdisk
_isdisk(::SparseMPO) = false
_isdisk(::RefMPO) = false
_isdisk(::RefMPS) = false

Base.size(t::DenseMPOTensor{4}) = map(dim, t.A |> x -> (codomain(x)[2], domain(x)[1]))
Base.size(::SparseMPOTensor{DL,D,DR,T}) where {DL,D,DR,T} = DL, DR

# length 定义一次, 覆盖全部容器 (Dense/Adjoint/Ref/Sparse 与 MPS/MPO 两侧)。
for W in (:DenseMPO, :AdjointMPO, :DenseMPS, :AdjointMPS, :RefMPS, :RefMPO, :SparseMPO)
    @eval Base.length(::($W){L}) where {L} = L
end

const _Indexable = Union{DenseMPO,AdjointMPO,DenseMPS,AdjointMPS,SparseMPO}

Base.firstindex(obj::_Indexable) = 1
Base.lastindex(obj::_Indexable) = lastindex(obj.ts)
Base.size(obj::_Indexable) = (lastindex(obj),)
Base.axes(obj::_Indexable) = Base.OneTo(lastindex(obj))

Base.firstindex(obj::RefMPS) = 1
Base.lastindex(obj::RefMPS) = lastindex(obj.ts)
Base.size(obj::RefMPS) = (lastindex(obj),)
Base.axes(obj::RefMPS) = Base.OneTo(lastindex(obj))

Base.firstindex(obj::RefMPO) = 1
Base.lastindex(obj::RefMPO) = lastindex(obj.ts)
Base.size(obj::RefMPO) = (lastindex(obj),)
Base.axes(obj::RefMPO) = Base.OneTo(lastindex(obj))

# 磁盘对象: getindex 反序列化 / setindex! 序列化; 内存对象: 引用直通。
Base.getindex(obj::_Indexable, i::Int64) = _isdisk(obj) ? (@timeit _local_io_timer() "deserialize" obj.ts[i]) : obj.ts[i]
Base.getindex(obj::_Indexable, stp::UnitRange) = _isdisk(obj) ? (@timeit _local_io_timer() "deserialize" [obj.ts[i] for i in stp]) : [obj.ts[i] for i in stp]
Base.getindex(obj::_Indexable, ::Colon) = _isdisk(obj) ? [obj.ts[i] for i in 1:length(obj.ts)] : obj.ts[:]

Base.setindex!(obj::_Indexable, val, i::Int64) = _isdisk(obj) ? (@timeit _local_io_timer() "serialize" obj.ts[i] = val) : (obj.ts[i] = val)
Base.setindex!(obj::_Indexable, vals, stp::UnitRange) = _isdisk(obj) ? (@timeit _local_io_timer() "serialize" for (i, v) in zip(stp, vals); obj.ts[i] = v; end) : (for (i, v) in zip(stp, vals); obj.ts[i] = v; end)
Base.setindex!(obj::_Indexable, vals, ::Colon) = _isdisk(obj) ? (for (i, v) in enumerate(vals); obj.ts[i] = v; end) : (obj.ts[:] = vals)

# Ref 层: mapping 逐元素作用 (默认 adjoint/identity); setindex! 空操作。
Base.getindex(obj::T, i::Int64) where {T<:Union{RefMPO,RefMPS}} = obj.mapping(obj.ts[i])
Base.getindex(obj::T, stp::UnitRange) where {T<:Union{RefMPO,RefMPS}} = obj.mapping.([obj.ts[i] for i in stp])
Base.getindex(obj::T, ::Colon) where {T<:Union{RefMPO,RefMPS}} = obj.mapping.(obj.ts[:])
Base.setindex!(::RefMPO, val, i::Int64) = nothing
Base.setindex!(::RefMPS, val, i::Int64) = nothing
Base.setindex!(::RefMPS, vals, stp::UnitRange) = nothing
Base.setindex!(::RefMPS, vals, ::Colon) = nothing

# ---------- MP 层容器: 范数 / 归一化 ----------
# norm 与 normalize! 都以 center 为单点 (Gauge 中心) 计算; normalize! 就地并返回 norm。

for W in (:DenseMPO, :DenseMPS, :AdjointMPO, :AdjointMPS)
    @eval function TensorKit.LinearAlgebra.normalize!(obj::($W){L}) where L
        @assert (site = obj.center[1]) == obj.center[2]
        t = obj[site]           # 磁盘对象: 反序列化; 内存对象: 引用
        n = normalize!(t)       # 就地归一化
        obj[site] = t           # 磁盘对象: 序列化写回; 内存对象: 无操作
        return n
    end
    @eval function TensorKit.norm(obj::($W){L}) where L
        @assert (site = obj.center[1]) == obj.center[2]
        return norm(obj[site])
    end
end

function TensorKit.LinearAlgebra.normalize!(obj::RefMPS)
    @assert (site = obj.center[1]) == obj.center[2]
    return normalize!(obj[site])  # RefMPS 的 setindex! 是空操作, 无需写回
end

function TensorKit.LinearAlgebra.normalize!(obj::RefMPO)
    @assert (site = obj.center[1]) == obj.center[2]
    return normalize!(obj[site])
end

function TensorKit.norm(obj::RefMPS)
    @assert (site = obj.center[1]) == obj.center[2]
    return norm(obj[site])
end

function TensorKit.norm(obj::RefMPO)
    @assert (site = obj.center[1]) == obj.center[2]
    return norm(obj[site])
end

# ---------- 稀疏环境张量: 索引 / 长度 / 迭代 ----------

Base.getindex(obj::T, i::Int64) where {T<:Union{SparseLeftEnvironmentTensor,SparseRightEnvironmentTensor}} = obj.A[i]
Base.getindex(obj::T, inds::AbstractVector{Int64}) where {T<:Union{SparseLeftEnvironmentTensor,SparseRightEnvironmentTensor}} = [obj.A[i] for i in inds]
Base.getindex(obj::T, i::Int...) where {T<:Union{SparseLeftEnvironmentTensor,SparseRightEnvironmentTensor}} = obj.A[i...]

Base.setindex!(obj::SparseLeftEnvironmentTensor, A::LeftEnvironmentTensor, i::Int64) = (obj.A[i] = A)
Base.setindex!(obj::SparseRightEnvironmentTensor, A::RightEnvironmentTensor, i::Int64) = (obj.A[i] = A)
Base.setindex!(obj::T, val, i::Int...) where {T<:Union{SparseLeftEnvironmentTensor,SparseRightEnvironmentTensor}} = (obj.A[i...] = val)

Base.length(obj::T) where {T<:Union{SparseLeftEnvironmentTensor,SparseRightEnvironmentTensor}} = length(obj.A)
Base.size(obj::T) where {T<:Union{SparseLeftEnvironmentTensor,SparseRightEnvironmentTensor}} = size(obj.A)
Base.iterate(obj::T, args...) where {T<:Union{SparseLeftEnvironmentTensor,SparseRightEnvironmentTensor}} = iterate(obj.A, args...)

# ---------- Tuple 混合向量的向量代数 (CBE/Lanczos 的 (tensor, scalar) 状态向量) ----------

TensorKit.scale(t::Tuple{T₁,T₂}) where {T₁<:AbstractTensorWrapper,T₂<:Number} = scale((t[1].A, t[2]))
TensorKit.scale!!(t::Tuple{T₁,T₂}) where {T₁<:AbstractTensorWrapper,T₂<:Number} = scale!!((t[1].A, t[2]))
TensorKit.scale!!(t::Tuple{T₁,T₂,T₃}) where {T₁<:AbstractTensorWrapper,T₂<:AbstractTensorWrapper,T₃<:Number} = scale!!((t[1].A, t[2].A, t[3]))
TensorKit.zerovector(t::Tuple{T₁,T₂}) where {T₁<:AbstractTensorWrapper,T₂<:Number} = zerovector((t[1].A, t[2]))
TensorKit.add!!(t::Tuple) = add!!(map(x -> x isa AbstractTensorWrapper ? x.A : x, t))
