
#基本思路，计算每一时刻的MPS，在update处插入演化与反向演化
#本质上来说不用求解本征值，直接对该点的MPS进行演化
function sweepTDVP1(ψ::Vector,H::Vector,
    t::Number,Nt::Int64,
    D_MPS::Int64)
    
    lsψ = Vector{Vector}(undef,Nt)
    lst = collect(range(0,t,Nt))
    τ = t/(Nt-1)/2

    lsψ[1] = ψ

    for i in 2:Nt
        println("evolution $i begin, t = $(round(lst[i];digits=3))")

        EnvR = RightEnv(ψ,H,1)
        ψ[1:2] .= RightUpdateTDVP1(ψ[1:2],H[1],EnvR,τ,D_MPS)
        EnvL = LeftEnv(ψ[1],H[1])

        for iL in 2:length(H)-1
            EnvR = RightEnv(ψ,H,iL)
            ψ[iL:iL+1] .= RightUpdateTDVP1(ψ[iL:iL+1],H[iL],EnvL,EnvR,τ,D_MPS)
            EnvL = PushRight(EnvL,ψ[iL],H[iL])
        end
        println("right sweep finished")

        EnvL = LeftEnv(ψ,H,1)
        ψ[end-1:end] .= LeftUpdateTDVP1(ψ[end-1:end],H[end],EnvL,2*τ,D_MPS)
        EnvR = RightEnv(ψ[end],H[end])

        for iL in length(H)-1:-1:2
            EnvL = LeftEnv(ψ,H,iL)
            ψ[iL-1:iL] .= LeftUpdateTDVP1(ψ[iL-1:iL],H[iL],EnvL,EnvR,τ,D_MPS)
            EnvR = PushLeft(EnvR,ψ[iL],H[iL])
        end
        println("left sweep finished")

        lsψ[i] = ψ
    end

    return lsψ,lst
end

#仿照DMRG，建立左右Move的函数

function RightUpdateTDVP1(ψs::Vector,Hi::AbstractTensorMap,
    EnvR::AbstractTensorMap,τ::Number,
    D_MPS::Int64)

    reduH1 = ReduHam1(Hi,EnvR)
    Aτ = Apply(ψs[1],EvolveOpr(reduH1,τ))

    Στ,thisMPS = RightSVD(Aτ,D_MPS)

    reduH0 = ReduHam1(LeftEnv(thisMPS,Hi),EnvR)
    Σ = Apply(Στ,EvolveOpr(reduH0,-τ))

    return RightMerge(Σ,thisMPS,ψs[2])
end

function RightUpdateTDVP1(ψs::Vector,Hi::AbstractTensorMap,
    EnvL::AbstractTensorMap,EnvR::AbstractTensorMap,τ::Number,
    D_MPS::Int64)

    reduH1 = ReduHam1(Hi,EnvL,EnvR)
    Aτ = Apply(ψs[1],EvolveOpr(reduH1,τ))

    Στ,thisMPS = RightSVD(Aτ,D_MPS)

    reduH0 = ReduHam1(PushRight(EnvL,thisMPS,Hi),EnvR)
    Σ = Apply(Στ,EvolveOpr(reduH0,-τ))

    return RightMerge(Σ,thisMPS,ψs[2])
end


function LeftUpdateTDVP1(ψs::Vector,Hi::AbstractTensorMap,
    EnvL::AbstractTensorMap,τ::Number,
    D_MPS::Int64)

    reduH1 = ReduHam1(Hi,EnvL)
    Aτ = Apply(ψs[2],EvolveOpr(reduH1,τ))

    Στ,thisMPS = LeftSVD(Aτ,D_MPS)

    reduH0 = ReduHam1(EnvL,RightEnv(thisMPS,Hi))
    Σ = Apply(Στ,EvolveOpr(reduH0,-τ))

    MPSs = LeftMerge(Σ,thisMPS,ψs[1])

    return MPSs
end

function LeftUpdateTDVP1(ψs::Vector,Hi::AbstractTensorMap,
    EnvL::AbstractTensorMap,EnvR::AbstractTensorMap,τ::Number,
    D_MPS::Int64)

    reduH1 = ReduHam1(Hi,EnvL,EnvR)
    Aτ = LeftApply(ψs[2],EvolveOpr(reduH1,τ))
    Aτ = permute(Aτ,(),(2,3,1))
    Στ,thisMPS = LeftSVD(Aτ,D_MPS)

    reduH0 = ReduHam1(EnvL,PushLeft(EnvR,thisMPS,Hi))
    Σ = Apply(Στ,EvolveOpr(reduH0,-τ))

    return LeftMerge(Σ,thisMPS,ψs[1])
end

function EvolveOpr(Ham::AbstractTensorMap,τ::Number)
    return exp(-1im * Ham * τ)
end

function Apply(ψ::AbstractTensorMap{ComplexSpace,1,2},Opr::AbstractTensorMap{ComplexSpace,2,4})
    @tensor ψτ[-1,-2,-3] ≔ ψ[1,2,3]*Opr[1,-1,2,3,-2,-3]
    return permute(ψτ,(1,),(2,3))
end

function Apply(ψ::AbstractTensorMap{ComplexSpace,1,1},Opr::AbstractTensorMap{ComplexSpace,2,2})
    @tensor ψτ[-1,-2] ≔ ψ[1,2]*Opr[1,-1,2,-2]
    return permute(ψτ,(1,),(2,))
end

function Apply(ψ::AbstractTensorMap{ComplexSpace,0,2},Opr::AbstractTensorMap{ComplexSpace,2,2})
    @tensor ψτ[-1,-2] ≔ ψ[1,2]*Opr[1,2,-1,-2]
    return permute(ψτ,(),(1,2))
end

function LeftApply(ψ::AbstractTensorMap{ComplexSpace,0,2},Opr::AbstractTensorMap{ComplexSpace,2,2})
    @tensor ψτ[-1,-2] ≔ ψ[1,2]*Opr[1,2,-1,-2]
    return permute(ψτ,(),(1,2))
end

function Apply(ψ::AbstractTensorMap{ComplexSpace,0,3},Opr::AbstractTensorMap{ComplexSpace,3,3})
    @tensor ψτ[-1,-2,-3] ≔ ψ[1,2,3]*Opr[1,2,3,-1,-2,-3]
    return permute(ψτ,(),(1,2,3))
end

function LeftApply(ψ::AbstractTensorMap{ComplexSpace,0,3},Opr::AbstractTensorMap{ComplexSpace,3,3})
    @tensor ψτ[-3,-1,-2] ≔ ψ[3,1,2]*Opr[1,2,3,-1,-2,-3]
    return permute(ψτ,(),(1,2,3))
end

function ReduHam0(EnvL::AbstractTensorMap{ComplexSpace,2,1},EnvR::AbstractTensorMap{ComplexSpace,1,2})
    @tensor reduH0[-1,-2,-3,-4] ≔ EnvL[-1,1,-3] * EnvR[-2,1,-4]
    return permute(reduH0,(1,2),(3,4))
end
