function RightEnv(ψ::Vector,H::Vector,site::Int64=1)

    EnvR = InitialRightEnv()
    for iL in length(ψ):-1:site+1
        EnvR = PushLeft(EnvR,ψ[iL],H[iL])
    end

    return EnvR
end

function RightLsEnv(ψ::Vector,H::Vector,site::Int64=1)

    LR = length(ψ) + 1 - site
    lsEnvR = Vector{AbstractTensorMap}(undef,LR)

    lsEnvR[LR] = InitialRightEnv()
    for i in LR-1:-1:1
        lsEnvR[i] = PushLeft(lsEnvR[i+1],ψ[site + i],H[site + i])
    end

    return lsEnvR
end

function LeftEnv(ψ::Vector,H::Vector,site::Int64=1)

    EnvL = InitialLeftEnv()
    for iL in 1:site-1
        EnvL = PushRight(EnvL,ψ[iL],H[iL])
    end

    return EnvL
end

function LeftLsEnv(ψ::Vector,H::Vector,site::Int64)

    LL = site
    lsEnvL = Vector{AbstractTensorMap}(undef,LL)

    lsEnvL[1] = InitialLeftEnv()
    for i in 2:LL
        lsEnvL[i] = PushRight(lsEnvL[i - 1],ψ[i - 1],H[i - 1])
    end

    return lsEnvL
end

function InitialRightEnv(space=ℂ)
    EnvR = TensorMap(reshape([1.0 + 0.0im],1,1,1), space^1, space^1 ⊗ space^1)
    return EnvR
end


function PushLeft(EnvR::AbstractTensorMap{ComplexSpace,1,2},ψi::AbstractTensorMap{ComplexSpace,1,2},Hi::AbstractTensorMap{ComplexSpace,2,2})
    @tensor EnvRR[-1,-2,-3] ≔ ψi[-1,4,1]*Hi[3,4,-2,5]*ψi'[5,2,-3]*EnvR[1,3,2]
    return permute(EnvRR,(1,),(2,3))
end

function InitialLeftEnv(space=ℂ)
    EnvL = TensorMap(reshape([1.0 + 0.0im],1,1,1), space^1 ⊗ space^1, space^1 )
    return EnvL
end

function PushRight(EnvL::AbstractTensorMap{ComplexSpace,2,1},ψi::AbstractTensorMap{ComplexSpace,1,2},Hi::AbstractTensorMap{ComplexSpace,2,2})
    @tensor EnvLL[-1,-2,-3] ≔ ψi[-1,1,4]*Hi[-2,4,3,5]*ψi'[2,5,-3]*EnvL[1,3,2]
    return permute(EnvLL,(1,2),(3,))
end
