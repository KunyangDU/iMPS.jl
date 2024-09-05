
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


function sweepDMRG1(ψ::Vector,H::Vector,Nsweep::Int64,LanczosLevel::Int64,D_MPS::Int64)
    lsE = Vector{Float64}(undef,Nsweep)

    # 部分改变单侧Env
    # 存储ψ依然使用全局存储

    for i in 1:Nsweep
        println("sweep $i begin")
        start_time = time()

        Eg = 0

        Eg, ψ[1:2] = RightUpdateDMRG1(ψ,H,1,LanczosLevel,D_MPS)
        EnvL = LeftMostEnv(ψ[1],H[1])
        for iL in 2:L-1
            EnvR = RightEnv(ψ,H,iL)
            reduH = InnerReduHam1(H[iL],EnvL,EnvR)

            Eg,Ev = groundEig(reduH,LanczosLevel)

            #ψ[iL:iL+1] = DirectedSVD(ψ,Ev',iL,"right",D_MPS)
            ψ[iL:iL+1] = RightMove(ψ[iL+1],Ev',D_MPS)

            EnvL = PushRight(EnvL,ψ[iL],H[iL])
        end
        println("right sweep finished")

        Eg, ψ[end-1:end] = LeftUpdateDMRG1(ψ,H,length(ψ),LanczosLevel,D_MPS)
        EnvR = RightMostEnv(ψ[end],H[end])
        for iL in L-1:-1:2
            EnvL = LeftEnv(ψ,H,iL)
            reduH = InnerReduHam1(H[iL],EnvL,EnvR)

            Eg,Ev = groundEig(reduH,LanczosLevel)

            #ψ[iL-1:iL] = DirectedSVD(ψ,Ev',iL,"left",D_MPS)
            ψ[iL-1:iL] = LeftMove(ψ[iL-1],Ev',D_MPS)

            EnvR = PushLeft(EnvR,ψ[iL],H[iL])
        end
        println("left sweep finished")
        println("sweep $i finished, Eg = $Eg, time consumed $(round(time()-start_time;digits=2))")
        lsE[i] = Eg
    end
    return ψ,lsE
end

function RightUpdateDMRG1(ψ::Vector,H::Vector,site::Int64,level::Int64,D_MPS::Int64)

    reduH = reduHam1(ψ,H,site)
    Eg,Ev = groundEig(reduH,level)
    MPSs = RightMove(ψ[site+1],Ev',site,D_MPS)

    return Eg, MPSs
end

function LeftUpdateDMRG1(ψ::Vector,H::Vector,site::Int64,level::Int64,D_MPS::Int64)

    reduH = reduHam1(ψ,H,site)
    Eg,Ev = groundEig(reduH,level)
    MPSs = LeftMove(ψ[site-1],Ev',site,D_MPS)

    return Eg, MPSs
end
