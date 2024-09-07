
#基本思路，计算每一时刻的MPS，在update处插入演化与反向演化
#本质上来说不用求解本征值，直接对该点的MPS进行演化

# 脚本拆分
# 注意Apply的冗余和EffHam的广泛性，可以简化代码

function sweepTDVP1(ψ::Vector,H::Vector,
    t::Number,Nt::Int64,
    D_MPS::Int64)
    
    lsψ = Vector{Vector}(undef,Nt)
    lst = collect(range(0,t,Nt))
    τ = t/(Nt-1)/2

    lsψ[1] = deepcopy(ψ)

    #MiMPO = TensorMap([1 0;0 -1],(ℂ^2)'→ (ℂ^2)')


    for i in 2:Nt
        #lsMi = Quant1(ψ,MiMPO,D_MPS)
        start_time = time()
        println("evolution $i begin, t = $(round(lst[i];digits=3))/J")

        EnvR = RightEnv(ψ,H,1)
        ψ[1:2] = RightUpdateTDVP1(ψ[1:2],H[1],EnvR,τ,D_MPS)
        EnvL = LeftEnv(ψ[1],H[1])

        for iL in 2:length(H)-1
            EnvR = RightEnv(ψ,H,iL)
            ψ[iL:iL+1] = RightUpdateTDVP1(ψ[iL:iL+1],H[iL],EnvL,EnvR,τ,D_MPS)
            EnvL = PushRight(EnvL,ψ[iL],H[iL])
        end
        println("right sweep finished")

        EnvL = LeftEnv(ψ,H,1)
        ψ[end-1:end] = LeftUpdateTDVP1(ψ[end-1:end],H[end],EnvL,2*τ,D_MPS)
        EnvR = RightEnv(ψ[end],H[end])

        for iL in length(H)-1:-1:2
            EnvL = LeftEnv(ψ,H,iL)
            ψ[iL-1:iL] = LeftUpdateTDVP1(ψ[iL-1:iL],H[iL],EnvL,EnvR,τ,D_MPS)
            EnvR = PushLeft(EnvR,ψ[iL],H[iL])
        end
        println("left sweep finished")

        println("evolution $i finished, time consumed $(round(time()-start_time;digits=2))s")

        lsψ[i] = deepcopy(ψ)
    end

    return lsψ,lst
end

#仿照DMRG，建立左右Move的函数

function RightUpdateTDVP1(ψs::Vector,Hi::AbstractTensorMap,
    EnvR::AbstractTensorMap,τ::Number,
    D_MPS::Int64)

    reduH1 = EffHam(Hi,EnvR)
    Aτ = Apply(ψs[1],EvolveOpr(reduH1,τ))

    Στ,thisMPS = RightSVD(Aτ,D_MPS)

    reduH0 = EffHam(LeftEnv(thisMPS,Hi),EnvR)
    Σ = Apply(Στ,EvolveOpr(reduH0,-τ))

    return collect(RightMerge(Σ,thisMPS,ψs[2]))
end

function RightUpdateTDVP1(ψs::Vector,Hi::AbstractTensorMap,
    EnvL::AbstractTensorMap,EnvR::AbstractTensorMap,τ::Number,
    D_MPS::Int64)

    reduH1 = EffHam(Hi,EnvL,EnvR)
    Aτ = Apply(ψs[1],EvolveOpr(reduH1,τ))

    Στ,thisMPS = RightSVD(Aτ,D_MPS)

    reduH0 = EffHam(PushRight(EnvL,thisMPS,Hi),EnvR)
    Σ = Apply(Στ,EvolveOpr(reduH0,-τ))

    return collect(RightMerge(Σ,thisMPS,ψs[2]))
end


function LeftUpdateTDVP1(ψs::Vector,Hi::AbstractTensorMap,
    EnvL::AbstractTensorMap,τ::Number,
    D_MPS::Int64)

    reduH1 = EffHam(Hi,EnvL)
    Aτ = Apply(ψs[2],EvolveOpr(reduH1,τ))

    Στ,thisMPS = LeftSVD(Aτ,D_MPS)

    reduH0 = EffHam(EnvL,RightEnv(thisMPS,Hi))
    Σ = Apply(Στ,EvolveOpr(reduH0,-τ))

    MPSs = collect(LeftMerge(Σ,thisMPS,ψs[1]))

    return MPSs
end

function LeftUpdateTDVP1(ψs::Vector,Hi::AbstractTensorMap,
    EnvL::AbstractTensorMap,EnvR::AbstractTensorMap,τ::Number,
    D_MPS::Int64)

    reduH1 = EffHam(Hi,EnvL,EnvR)
    Aτ = LeftApply(ψs[2],EvolveOpr(reduH1,τ))
    Aτ = permute(Aτ,(),(2,3,1))
    Στ,thisMPS = LeftSVD(Aτ,D_MPS)

    reduH0 = EffHam(EnvL,PushLeft(EnvR,thisMPS,Hi))
    Σ = Apply(Στ,EvolveOpr(reduH0,-τ))

    return collect(LeftMerge(Σ,thisMPS,ψs[1]))
end

