
#################### INITIAL ENVIRONMENT ####################

function InitialEnv(lsEnv::Vector{AbstractTensorMap{ComplexSpace,codN,dN}};
    space=ℂ, ds::Vector{Int64} = ones(Int64,codN + dN)) where {codN,dN}
    return TensorMap(ones(ComplexF64,prod(ds)), ⊗([space^(ds[i]) for i in 1:codN]...) , ⊗([space^(ds[codN + i]) for i in 1:dN]...) )
end

function InitialEnv(lsEnv::Vector{AbstractTensorMap{ComplexSpace,0,dN}};
    space=ℂ, ds::Vector{Int64} = ones(Int64,dN)) where dN
    return Tensor(ones(ComplexF64,prod(ds)) , ⊗([space^(ds[i]) for i in 1:dN]...) ) |> x -> permute(x,(),Tuple(1:dN))
end

function InitialEnv(lsEnv::Vector{AbstractTensorMap{ComplexSpace,codN,0}};
    space=ℂ, ds::Vector{Int64} = ones(Int64,codN)) where codN
    return Tensor(ones(ComplexF64,prod(ds)), ⊗([space^(ds[i]) for i in 1:codN]...)) |> x -> permute(x,Tuple(1:codN),())
end


function InitialRightEnv(;space=ℂ,order::Int64=3)
    if order == 3
        EnvR = TensorMap(reshape([1.0 + 0.0im],1,1,1), space^1, space^1 ⊗ space^1)
    elseif order == 2
        EnvR = TensorMap(reshape([1.0 + 0.0im],1,1), space^1, space^1)
    end
    return EnvR
end

function InitialRightEnv(lsEnvR::Vector{AbstractTensorMap{ComplexSpace,codN,dN}};space=ℂ) where {codN,dN}
    return TensorMap(reshape(ones(ComplexF64,codN*dN),codN,dN), ⊗(space^1 for _ in 1:codN) , ⊗(space^1 for _ in 1:dN) )
end


function InitialLeftEnv(;space=ℂ,order::Int64=3)
    if order == 3
        EnvL = TensorMap(reshape([1.0 + 0.0im],1,1,1), space^1 ⊗ space^1, space^1 )
    elseif order == 2
        EnvL = TensorMap(reshape([1.0 + 0.0im],1,1), space^1, space^1)
    end

    return EnvL
end

function InitialLeftEnv(lsEnvL::Vector{AbstractTensorMap{ComplexSpace,codN,dN}};space=ℂ) where {codN,dN}
    return TensorMap(reshape(ones(ComplexF64,codN*dN),codN,dN), ⊗([space^1 for i in 1:codN]...) , ⊗([space^1 for i in 1:dN]...) )
end

