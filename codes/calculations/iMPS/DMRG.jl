
function sweepDMRG1(ψ::Vector,H::Vector,Nsweep::Int64,LanczosLevel::Int64,D_MPS::Int64)
    lsE = Vector{Float64}(undef,Nsweep)

    # 部分改变单侧Env
    # 存储ψ依然使用全局存储
    # 考虑端点处补充虚拟张量？

    for i in 1:Nsweep
        println("sweep $i begin")
        start_time = time()

        Eg = 0

        EnvR = RightEnv(ψ,H,1)
        Eg, ψ[1:2] = RightUpdateDMRG1(ψ[2],H[1],EnvR,LanczosLevel,D_MPS)
        EnvL = LeftEnv(ψ[1],H[1])
        for iL in 2:L-1
            EnvR = RightEnv(ψ,H,iL)
            Eg,ψ[iL:iL+1] = RightUpdateDMRG1(ψ[iL+1],H[iL],EnvL,EnvR,LanczosLevel,D_MPS)
            EnvL = PushRight(EnvL,ψ[iL],H[iL])
        end
        println("right sweep finished")

        EnvL = LeftEnv(ψ,H,1)
        Eg, ψ[end-1:end] = LeftUpdateDMRG1(ψ[end-1],H[end],EnvL,LanczosLevel,D_MPS)
        EnvR = RightEnv(ψ[end],H[end])
        for iL in L-1:-1:2
            EnvL = LeftEnv(ψ,H,iL)
            Eg,ψ[iL-1:iL] = LeftUpdateDMRG1(ψ[iL-1],H[iL],EnvL,EnvR,LanczosLevel,D_MPS)
            EnvR = PushLeft(EnvR,ψ[iL],H[iL])
        end
        println("left sweep finished")
        println("sweep $i finished, Eg = $Eg, time consumed $(round(time()-start_time;digits=2))")
        lsE[i] = Eg
    end
    return ψ,lsE
end

function RightUpdateDMRG1(nextψ::AbstractTensorMap,Hi::AbstractTensorMap,
    EnvR::AbstractTensorMap,
    LanczosLevel::Int64,D_MPS::Int64)

    reduH = ReduHam1(Hi,EnvR)
    Eg,Ev = groundEig(reduH,LanczosLevel)
    MPSs = RightMove(nextψ,Ev',D_MPS)

    return Eg, MPSs
end

function RightUpdateDMRG1(nextψ::AbstractTensorMap,Hi::AbstractTensorMap,
    EnvL::AbstractTensorMap,EnvR::AbstractTensorMap,
    LanczosLevel::Int64,D_MPS::Int64)

    reduH = ReduHam1(Hi,EnvL,EnvR)
    Eg,Ev = groundEig(reduH,LanczosLevel)
    MPSs = RightMove(nextψ,Ev',D_MPS)

    return Eg, MPSs
end

function LeftUpdateDMRG1(nextψ::AbstractTensorMap,Hi::AbstractTensorMap,
    EnvL::AbstractTensorMap,EnvR::AbstractTensorMap,
    LanczosLevel::Int64,D_MPS::Int64)

    reduH = ReduHam1(Hi,EnvL,EnvR)
    Eg,Ev = groundEig(reduH,LanczosLevel)
    MPSs = LeftMove(nextψ,Ev',D_MPS)

    return Eg, MPSs
end

function LeftUpdateDMRG1(nextψ::AbstractTensorMap,Hi::AbstractTensorMap,
    EnvL::AbstractTensorMap,
    LanczosLevel::Int64,D_MPS::Int64)

    reduH = ReduHam1(Hi,EnvL)
    Eg,Ev = groundEig(reduH,LanczosLevel)
    MPSs = LeftMove(nextψ,Ev',D_MPS)

    return Eg, MPSs
end

# stable

function RightUpdateDMRG1(ψ::Vector,H::Vector,site::Int64,LanczosLevel::Int64,D_MPS::Int64)

    reduH = ReduHam1(ψ,H,site)
    Eg,Ev = groundEig(reduH,LanczosLevel)
    MPSs = RightMove(ψ[site+1],Ev',D_MPS)

    return Eg, MPSs
end

function LeftUpdateDMRG1(ψ::Vector,H::Vector,site::Int64,LanczosLevel::Int64,D_MPS::Int64)

    reduH = ReduHam1(ψ,H,site)
    Eg,Ev = groundEig(reduH,LanczosLevel)
    MPSs = LeftMove(ψ[site-1],Ev',D_MPS)

    return Eg, MPSs
end


#= function sweepDMRG1(ψ::Vector,H::Vector,Nsweep::Int64,LanczosLevel::Int64,D_MPS::Int64)
    lsE = Vector{Float64}(undef,Nsweep)

    for i in 1:Nsweep
        println("sweep $i begin")
        start_time = time()
        Eg = 0
        for iL in 1:L-1
            showdomain(ψ[iL+1])
            Eg, ψ[iL:iL+1] = updateDMRG1(ψ,H,iL,LanczosLevel,"right",D_MPS)
            showdomain(ψ[iL])
        end
        println("right sweep finished")
        for iL in L:-1:2
            Eg, ψ[iL-1:iL] = updateDMRG1(ψ,H,iL,LanczosLevel,"left",D_MPS)
        end
        println("left sweep finished")
        println("sweep $i finished, Eg = $Eg, time consumed $(round(time()-start_time;digits=2))")
        lsE[i] = Eg
    end
    return ψ,lsE
end =#
