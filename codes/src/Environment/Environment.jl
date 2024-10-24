
#################### UNILATERAL ENVIRONMENT ####################

function RightEnv(ψ::Vector,H::Vector{AbstractTensorMap{ComplexSpace,2,2}},site::Int64)

    EnvR = InitialEnv(Vector{AbstractTensorMap{ComplexSpace,1,2}}(undef,1);ds = [dim(domain(ψ[end])[2]), dim(codomain(H[end])[1]),dim(domain(ψ[end])[2])])
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

    EnvL = InitialEnv(Vector{AbstractTensorMap{ComplexSpace,2,1}}(undef,1);ds = [dim(domain(ψ[1])[1]), dim(domain(H[1])[1]),dim(domain(ψ[1])[1])])
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
    ψ::Vector{Union{AbstractTensorMap{ComplexSpace,1,2},AbstractTensorMap{ComplexSpace,0,3}}},
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

function RightLsEnv(
    ψ1::Vector{Union{AbstractTensorMap{ComplexSpace,1,2},AbstractTensorMap{ComplexSpace,0,3}}},
    H::Vector{AbstractTensorMap{ComplexSpace,2,2}},
    ψ2::Vector{Union{AbstractTensorMap{ComplexSpace,1,2},AbstractTensorMap{ComplexSpace,0,3}}}
    ,site::Int64)

    LR = length(H) + 1 - site
    lsEnvR = Vector{AbstractTensorMap}(undef,LR)

    lsEnvR[LR] = InitialRightEnv()
    
    for i in LR-1:-1:1
        lsEnvR[i] = PushLeft(lsEnvR[i+1],ψ1[site + i],H[site + i],ψ2[site + i])
    end

    return lsEnvR
end


function LeftLsEnv(
    ψ::Vector{Union{AbstractTensorMap{ComplexSpace,1,2},AbstractTensorMap{ComplexSpace,0,3}}},
    H::Vector{AbstractTensorMap{ComplexSpace,2,2}},
    site::Int64
    )

    LL = site
    lsEnvL = Vector{AbstractTensorMap}(undef,LL)

    lsEnvL[1] = InitialLeftEnv()
    for i in 2:LL
        lsEnvL[i] = PushRight(lsEnvL[i - 1],ψ[i - 1],H[i - 1])
    end

    return lsEnvL
end

function LeftLsEnv(
    ψ1::Vector{Union{AbstractTensorMap{ComplexSpace,1,2},AbstractTensorMap{ComplexSpace,0,3}}},
    H::Vector{AbstractTensorMap{ComplexSpace,2,2}},
    ψ2::Vector{Union{AbstractTensorMap{ComplexSpace,1,2},AbstractTensorMap{ComplexSpace,0,3}}},
    site::Int64
    )

    LL = site
    lsEnvL = Vector{AbstractTensorMap}(undef,LL)

    lsEnvL[1] = InitialLeftEnv()
    for i in 2:LL
        lsEnvL[i] = PushRight(lsEnvL[i - 1],ψ1[i - 1],H[i - 1],ψ2[i - 1])
    end

    return lsEnvL
end

############### normlized Opr ###############

function RightLsEnv(
    ψ::Vector{Union{AbstractTensorMap{ComplexSpace,0,3},AbstractTensorMap{ComplexSpace,1,2}}},
    Opr::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}},
    site::Int64
    )

    LR = length(ψ) + 1 - site
    lsEnvR = Vector{AbstractTensorMap{ComplexSpace,2,1}}(undef,LR)

    lsEnvR[LR] = InitialEnv(lsEnvR;ds = map(x -> dims(domain(x[end]))[2],[ψ,Opr,ψ]))

    for i in LR-1:-1:1
        lsEnvR[i] = PushLeft(lsEnvR[i+1],ψ[site + i],Opr[site + i])
    end

    return lsEnvR
end

function RightLsEnv(
    ψ1::Vector{Union{AbstractTensorMap{ComplexSpace,1,2},AbstractTensorMap{ComplexSpace,0,3}}},
    Opr::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}},
    ψ2::Vector{Union{AbstractTensorMap{ComplexSpace,1,2},AbstractTensorMap{ComplexSpace,0,3}}}
    ,site::Int64)

    LR = length(Opr) + 1 - site
    lsEnvR = Vector{AbstractTensorMap{ComplexSpace,2,1}}(undef,LR)

    lsEnvR[LR] = InitialEnv(lsEnvR;ds = map(x -> dims(domain(x[end]))[2],[ψ1,Opr,ψ2]))
    
    for i in LR-1:-1:1
        lsEnvR[i] = PushLeft(lsEnvR[i+1],ψ1[site + i],Opr[site + i],ψ2[site + i])
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

function LeftLsEnv(
    ψ::Vector{Union{AbstractTensorMap{ComplexSpace,0,3},AbstractTensorMap{ComplexSpace,1,2}}},
    Opr::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}},
    site::Int64
    )

    LL = site
    lsEnvL = Vector{AbstractTensorMap{ComplexSpace,2,1}}(undef,LL)

    lsEnvL[1] = InitialEnv(lsEnvL;ds = map(x -> dims(domain(x[1]))[1],[ψ,Opr,ψ]))
    
    for i in 2:LL
        lsEnvL[i] = PushRight(lsEnvL[i-1],ψ[site - i],Opr[site - i])
    end

    return lsEnvL
end

function LeftLsEnv(
    ψ1::Vector{Union{AbstractTensorMap{ComplexSpace,1,2},AbstractTensorMap{ComplexSpace,0,3}}},
    Opr::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}},
    ψ2::Vector{Union{AbstractTensorMap{ComplexSpace,1,2},AbstractTensorMap{ComplexSpace,0,3}}},
    site::Int64
    )

    LL = site
    lsEnvL = Vector{AbstractTensorMap{ComplexSpace,2,1}}(undef,LL)

    lsEnvL[1] = InitialEnv(lsEnvL;ds = map(x -> dims(domain(x[1]))[1],[ψ1,Opr,ψ2]))
    for i in 2:LL
        lsEnvL[i] = PushRight(lsEnvL[i - 1],ψ1[i - 1],Opr[i - 1],ψ2[i - 1])
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

