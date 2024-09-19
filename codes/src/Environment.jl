
#################### UNILATERAL ENVIRONMENT ####################

function RightEnv(ψ::Vector,H::Vector{AbstractTensorMap{ComplexSpace,2,2}},site::Int64)

    EnvR = InitialRightEnv()
    for iL in length(ψ):-1:site+1
        EnvR = PushLeft(EnvR,ψ[iL],H[iL])
    end

    return EnvR
end

function RightEnv(ψ1::Vector,H::Vector,ψ2::Vector,site::Int64)

    EnvR = InitialRightEnv()
    for iL in length(H):-1:site+1
        EnvR = PushLeft(EnvR,ψ1[iL],H[iL],ψ2[iL])
    end

    return EnvR
end

function RightEnv(ψ1::Vector,ψ2::Vector,site::Int64)

    EnvR = InitialRightEnv(;order=2)
    for iL in length(ψ2):-1:site+1
        EnvR = PushLeft(EnvR,ψ1[iL],ψ2[iL])
    end

    return EnvR
end

function RightEnv(
    ψ::Vector,
    H::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}},
    site::Int64
    )

    EnvR = InitialEnv(Vector{AbstractTensorMap{ComplexSpace,2,1}}(undef,1))
    for iL in length(ψ):-1:site+1
        EnvR = PushLeft(EnvR,ψ[iL],H[iL])
    end

    return EnvR
end

function RightEnv(
    Opr::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}},
    site::Int64
    )

    EnvR = InitialEnv(Vector{AbstractTensorMap{ComplexSpace,1,0}}(undef,1);ds = [dims(domain(Opr[end]))[2]])
    for iL in length(Opr):-1:site+1
        EnvR = PushLeft(EnvR,Opr[iL])
    end

    return EnvR
end

function RightEnv(
    Opr1::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}},
    Opr2::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}},
    site::Int64
    )

    EnvR = InitialEnv(Vector{AbstractTensorMap{ComplexSpace,1,1}}(undef,1);ds = map(x -> dims(domain(x[end]))[2],[Opr1,Opr2]))
    for iL in length(Opr1):-1:site+1
        EnvR = PushLeft(EnvR,Opr1[iL],Opr2[iL])
    end

    return EnvR
end

function RightEnv(
    Opr1::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}},
    Opr2::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}},
    Opr3::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}},
    site::Int64
    )

    EnvR = InitialEnv(Vector{AbstractTensorMap{ComplexSpace,2,1}}(undef,1);ds = map(x -> dims(domain(x[end]))[2],[Opr1,Opr2,Opr3]))
    for iL in length(Opr1):-1:site+1
        EnvR = PushLeft(EnvR,Opr1[iL],Opr2[iL],Opr3[iL])
    end

    return EnvR
end

function LeftEnv(
    ψ::Vector,
    H::Vector{AbstractTensorMap{ComplexSpace,2,2}},
    site::Int64
    )

    EnvL = InitialLeftEnv()
    for iL in 1:site-1
        EnvL = PushRight(EnvL,ψ[iL],H[iL])
    end

    return EnvL
end

function LeftEnv(
    ψ::Vector,
    H::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}},
    site::Int64
    )

    EnvL = Vector{AbstractTensorMap{ComplexSpace,2,1}}(undef,1)

    EnvL[1] = InitialEnv(EnvL)
    for iL in 1:site-1
        EnvL[1] = PushRight(EnvL,ψ[iL],H[iL])
    end

    return EnvL[1]
end

function LeftEnv(
    Opr::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}},
    site::Int64
    )

    EnvL = InitialEnv(Vector{AbstractTensorMap{ComplexSpace,1,0}}(undef,1);ds = [dims(domain(Opr[end]))[2]])
    for iL in 1:site-1
        EnvL = PushRight(EnvL,Opr[iL])
    end

    return EnvL
end

function LeftEnv(
    Opr1::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}},
    Opr2::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}},
    site::Int64
    )

    EnvL = InitialEnv(Vector{AbstractTensorMap{ComplexSpace,1,1}}(undef,1);ds = map(x -> dims(domain(x[end]))[2],[Opr1,Opr2]))
    for iL in 1:site-1
        EnvL = PushRight(EnvL,Opr1[iL],Opr2[iL])
    end

    return EnvL
end

function LeftEnv(
    Opr1::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}},
    Opr2::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}},
    Opr3::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}},
    site::Int64
    )

    EnvL = InitialEnv(Vector{AbstractTensorMap{ComplexSpace,2,1}}(undef,1);ds = map(x -> dims(domain(x[end]))[2],[Opr1,Opr2,Opr3]))
    for iL in 1:site-1
        EnvL = PushRight(EnvL,Opr1[iL],Opr2[iL],Opr3[iL])
    end

    return EnvL
end

function LeftEnv(ψ1::Vector,H::Vector,ψ2::Vector,site::Int64)

    EnvL = InitialLeftEnv()
    for iL in 1:site-1
        EnvL = PushRight(EnvL,ψ1[iL],H[iL],ψ2[iL])
    end

    return EnvL
end


#################### UNILATERAL ENVIRONMENT LIST ####################

function RightLsEnv(
    ψ::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}},
    H::Vector{AbstractTensorMap{ComplexSpace,2,2}},
    site::Int64
    )

    LR = length(ψ) + 1 - site
    lsEnvR = Vector{AbstractTensorMap}(undef,LR)

    lsEnvR[LR] = InitialRightEnv()
    for i in LR-1:-1:1
        lsEnvR[i] = PushLeft(lsEnvR[i+1],ψ[site + i],H[site + i])
    end

    return lsEnvR
end

function RightLsEnv(ψ1::Vector,H::Vector,ψ2::Vector,site::Int64)

    LR = length(H) + 1 - site
    lsEnvR = Vector{AbstractTensorMap}(undef,LR)

    lsEnvR[LR] = InitialRightEnv()
    
    for i in LR-1:-1:1
        lsEnvR[i] = PushLeft(lsEnvR[i+1],ψ1[site + i],H[site + i],ψ2[site + i])
    end

    return lsEnvR
end

function RightLsEnv(
    Opr1::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}},
    Opr2::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}},
    site::Int64
    )

    LR = length(Opr1) + 1 - site
    lsEnvR = Vector{AbstractTensorMap{ComplexSpace,1,1}}(undef,LR)

    lsEnvR[LR] = InitialEnv(lsEnvR;ds = map(x -> dims(domain(x[end]))[2],[Opr1,Opr2]))
    
    for i in LR-1:-1:1
        lsEnvR[i] = PushLeft(lsEnvR[i+1],Opr1[site + i],Opr2[site + i])
    end

    return lsEnvR
end

function RightLsEnv(
    Opr1::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}},
    Opr2::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}},
    Opr3::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}},
    site::Int64
    )

    LR = length(Opr1) + 1 - site
    lsEnvR = Vector{AbstractTensorMap{ComplexSpace,2,1}}(undef,LR)

    lsEnvR[LR] = InitialEnv(lsEnvR;ds = map(x -> dims(domain(x[end]))[2],[Opr1,Opr2,Opr3]))
    
    for i in LR-1:-1:1
        lsEnvR[i] = PushLeft(lsEnvR[i+1],Opr1[site + i],Opr2[site + i],Opr3[site + i])
    end

    return lsEnvR
end

function RightLsEnv(
    ψ1::Vector{Union{AbstractTensorMap{ComplexSpace,0,3},AbstractTensorMap{ComplexSpace,1,2}}},
    ψ2::Vector{Union{AbstractTensorMap{ComplexSpace,0,3},AbstractTensorMap{ComplexSpace,1,2}}},
    site::Int64
    )

    LR = length(ψ1) + 1 - site
    lsEnvR = Vector{AbstractTensorMap{ComplexSpace,1,1}}(undef,LR)

    lsEnvR[LR] = InitialEnv(lsEnvR;ds = map(x -> dims(domain(x[end]))[2],[ψ1,ψ2]))
    for i in LR-1:-1:1
        lsEnvR[i] = PushLeft(lsEnvR[i+1],ψ1[site + i],ψ2[site + i])
    end

    return lsEnvR
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

function LeftLsEnv(ψ1::Vector,H::Vector,ψ2::Vector,site::Int64)

    LL = site
    lsEnvL = Vector{AbstractTensorMap}(undef,LL)

    lsEnvL[1] = InitialLeftEnv()
    for i in 2:LL
        lsEnvL[i] = PushRight(lsEnvL[i - 1],ψ1[i - 1],H[i - 1],ψ2[i - 1])
    end

    return lsEnvL
end

function LeftLsEnv(
    Opr1::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}},
    Opr2::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}},
    site::Int64
    )

    LL = site
    lsEnvL = Vector{AbstractTensorMap{ComplexSpace,1,1}}(undef,LL)

    lsEnvL[1] = InitialEnv(lsEnvL;ds = map(x -> dims(domain(x[end]))[2],[Opr1,Opr2]))
    
    for i in 2:LL
        lsEnvL[i] = PushRight(lsEnvL[i-1],Opr1[site - i],Opr2[site - i])
    end

    return lsEnvL
end

function LeftLsEnv(
    Opr1::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}},
    Opr2::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}},
    Opr3::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}},
    site::Int64
    )

    LL = site
    lsEnvL = Vector{AbstractTensorMap{ComplexSpace,2,1}}(undef,LL)

    lsEnvL[1] = InitialEnv(lsEnvL;ds = map(x -> dims(domain(x[end]))[2],[Opr1,Opr2,Opr3]))
    
    for i in 2:LL
        lsEnvL[i] = PushRight(lsEnvL[i-1],Opr1[site - i],Opr2[site - i],Opr3[site - i])
    end

    return lsEnvL
end

function LeftLsEnv(
    ψ1::Vector{Union{AbstractTensorMap{ComplexSpace,1,2},AbstractTensorMap{ComplexSpace,0,3}}},
    ψ2::Vector{Union{AbstractTensorMap{ComplexSpace,1,2},AbstractTensorMap{ComplexSpace,0,3}}},
    site::Int64
    )

    LL = site
    lsEnvL = Vector{AbstractTensorMap{ComplexSpace,1,1}}(undef,LL)

    lsEnvL[1] = InitialEnv(lsEnvL;ds = map(x -> dims(domain(x[end]))[2],[ψ1,ψ2]))
    
    for i in 2:LL
        lsEnvL[i] = PushRight(lsEnvL[i-1],ψ1[site - i],ψ2[site - i])
    end

    return lsEnvL
end

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

#################### PUSH ####################

function PushLeft(EnvR::AbstractTensorMap{ComplexSpace,1,2},ψi::AbstractTensorMap{ComplexSpace,1,2},Hi::AbstractTensorMap{ComplexSpace,2,2})
    @tensor EnvRR[-1,-2,-3] ≔ ψi[-1,4,1]*Hi[3,4,-2,5]*ψi'[5,2,-3]*EnvR[1,3,2]
    return permute(EnvRR,(1,),(2,3))
end

function PushLeft(EnvR::AbstractTensorMap{ComplexSpace,1,2},ψi1::AbstractTensorMap{ComplexSpace,1,2},Hi::AbstractTensorMap{ComplexSpace,2,2},ψi2::AbstractTensorMap{ComplexSpace,1,2})
    @tensor EnvRR[-1,-2,-3] ≔ ψi1[-1,4,1]*Hi[3,4,-2,5]*ψi2'[5,2,-3]*EnvR[1,3,2]
    return permute(EnvRR,(1,),(2,3))
end

function PushLeft(EnvR::AbstractTensorMap{ComplexSpace,1,1},ψi1::AbstractTensorMap{ComplexSpace,1,2},ψi2::AbstractTensorMap{ComplexSpace,1,2})
    @tensor EnvRR[-1,-2] ≔ ψi1[-1,3,1] * ψi2'[3,2,-2] * EnvR[1,2]
    return permute(EnvRR,(1,),(2,))
end

function PushLeft(EnvR::AbstractTensorMap{ComplexSpace,2,1},ψi::AbstractTensorMap{ComplexSpace,1,2},Opri::AbstractTensorMap{ComplexSpace,2,2})
    @tensor EnvRR[-1,-2,-3] ≔ ψi[-1,3,1] * Opri[3,-2,5,2] * ψi'[5,4,-3] * EnvR[1,2,4]
    return permute(EnvRR,(1,2),(3,))
end

function PushLeft(EnvR::AbstractTensorMap{ComplexSpace,1,1},Opri1::AbstractTensorMap{ComplexSpace,2,2},Opri2::AbstractTensorMap{ComplexSpace,2,2})
    @tensor EnvRR[-1,-2] ≔ Opri1[4,-1,3,1] * Opri2'[3,2,4,-2] * EnvR[1,2]
    return permute(EnvRR,(1,),(2,))
end

function PushLeft(EnvR::AbstractTensorMap{ComplexSpace,2,1},Opri1::AbstractTensorMap{ComplexSpace,2,2},Opri2::AbstractTensorMap{ComplexSpace,2,2},Opri3::AbstractTensorMap{ComplexSpace,2,2})
    @tensor EnvRR[-1,-2,-3] ≔ Opri1[5,-1,3,1] * Opri2[3,-2,6,2] * Opri3'[6,4,5,-3] * EnvR[1,2,4]
    return permute(EnvRR,(1,2),(3,))
end

function PushLeft(EnvR::AbstractTensorMap{ComplexSpace,1,0},Opri::AbstractTensorMap{ComplexSpace,2,2})
    @tensor EnvRR[-1] ≔ Opri[2,-1,2,1] * EnvR[1]
    return permute(EnvRR,(1,),())
end

function PushRight(EnvL::AbstractTensorMap{ComplexSpace,2,1},ψi::AbstractTensorMap{ComplexSpace,1,2},Hi::AbstractTensorMap{ComplexSpace,2,2})
    @tensor EnvLL[-1,-2,-3] ≔ ψi[-1,1,4]*Hi[-2,4,3,5]*ψi'[2,5,-3]*EnvL[1,3,2]
    return permute(EnvLL,(1,2),(3,))
end

function PushRight(EnvL::AbstractTensorMap{ComplexSpace,2,1},ψi1::AbstractTensorMap{ComplexSpace,1,2},Hi::AbstractTensorMap{ComplexSpace,2,2},ψi2::AbstractTensorMap{ComplexSpace,1,2})
    @tensor EnvLL[-1,-2,-3] ≔ ψi1[-1,1,4]*Hi[-2,4,3,5]*ψi2'[2,5,-3]*EnvL[1,3,2]
    return permute(EnvLL,(1,2),(3,))
end

function PushRight(EnvL::AbstractTensorMap{ComplexSpace,1,1},Opri1::AbstractTensorMap{ComplexSpace,2,2},Opri2::AbstractTensorMap{ComplexSpace,2,2})
    @tensor EnvLL[-1,-2] ≔ Opri1[-1,3,1,4] * Opri2'[2,4,-2,3] * EnvL[1,2]
    return permute(EnvLL,(1,),(2,))
end

function PushRight(EnvL::AbstractTensorMap{ComplexSpace,1,1},ψi1::AbstractTensorMap{ComplexSpace,1,2},ψi2::AbstractTensorMap{ComplexSpace,1,2})
    @tensor EnvLL[-1,-2] ≔ ψi1[-1,1,3] * ψi2'[2,3,-2] * EnvL[1,2]
    return permute(EnvLL,(1,),(2,))
end

function PushRight(EnvL::AbstractTensorMap{ComplexSpace,2,1},Opri1::AbstractTensorMap{ComplexSpace,2,2},Opri2::AbstractTensorMap{ComplexSpace,2,2},Opri3::AbstractTensorMap{ComplexSpace,2,2})
    @tensor EnvLL[-1,-2,-3] ≔ Opri1[-1,5,1,3] * Opri2[-2,3,2,6] * Opri3'[4,6,-3,5] * EnvL[1,2,4]
    return permute(EnvLL,(1,2),(3,))
end

function PushRight(EnvL::AbstractTensorMap{ComplexSpace,1,0},Opri::AbstractTensorMap{ComplexSpace,2,2})
    @tensor EnvLL[-1] ≔ EnvL[1] * Opri[-1,2,1,2]
    return permute(EnvLL,(1,),())
end
