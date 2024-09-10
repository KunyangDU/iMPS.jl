
function sweepDMRG1(ψ::Vector,H::Vector,Nsweep::Int64,LanczosLevel::Int64,D_MPS::Int64)
    # 存储ψ依然使用全局存储

    L = length(H)

    lsE = Vector{Float64}(undef,Nsweep)

    # calculate the Environment
    lsEnv = vcat(LeftLsEnv(ψ,H,1),RightLsEnv(ψ,H,1))
    totaltruncerror = 0
    for i in 1:Nsweep
        temptruncerr = 0
        println("sweep $i")
        start_time = time()

        Eg = 0
        println(">>>>>> begin >>>>>>")
        for i in 1:L-1
            Eg,Ev = LocalEigen(H[i],lsEnv[i],lsEnv[i+1],LanczosLevel)
            ψ[i:i+1],temptruncerr = RightMove(ψ[i+1],Ev,D_MPS)
            #Eg, ψ[i:i+1] = RightUpdateDMRG1(ψ[i+1],H[i],lsEnv[i],lsEnv[i+1],LanczosLevel,D_MPS)
            lsEnv[i+1] = PushRight(lsEnv[i],ψ[i],H[i])
            totaltruncerror = max(totaltruncerror,temptruncerr)
        end
        println(">>>>>> finished >>>>>>")

        println("<<<<<< begin <<<<<<")
        for i in L:-1:2
            Eg,Ev = LocalEigen(H[i],lsEnv[i],lsEnv[i+1],LanczosLevel)
            ψ[i-1:i],temptruncerr = LeftMove(ψ[i-1],Ev,D_MPS)
            #Eg, ψ[i-1:i] = LeftUpdateDMRG1(ψ[i-1],H[i],lsEnv[i],lsEnv[i+1],LanczosLevel,D_MPS)
            lsEnv[i] = PushLeft(lsEnv[i+1],ψ[i],H[i])
            totaltruncerror = max(totaltruncerror,temptruncerr)
        end
        println("<<<<<< finished <<<<<<")

        println("sweep $i finished, Eg = $Eg, time consumed $(round(time()-start_time;digits=2)), max truncation error = $(totaltruncerror)")
        lsE[i] = Eg
    end

    return ψ,lsE
end


function sweepDMRG2(ψ::Vector,H::Vector,
    Nsweep::Int64,LanczosLevel::Int64,D_MPS::Int64)
    # 存储ψ依然使用全局存储

    L = length(H)

    lsE = Vector{Float64}(undef,Nsweep)

    # calculate the Environment
    lsEnv = vcat(LeftLsEnv(ψ,H,1),RightLsEnv(ψ,H,1))

    totaltruncerror = 0
    for i in 1:Nsweep
        temptruncerr = 0
        println("sweep $i")
        start_time = time()

        Eg = 0
        println(">>>>>> begin >>>>>>")
        for i in 1:L-1
            Eg,Ev = LocalEigen(H[i:i+1],lsEnv[i],lsEnv[i+2],LanczosLevel)
            ψ[i:i+1],temptruncerr = RightSVD(Ev,D_MPS)
            #Eg, ψ[i:i+1] = RightUpdateDMRG1(ψ[i+1],H[i],lsEnv[i],lsEnv[i+1],LanczosLevel,D_MPS)
            lsEnv[i+1] = PushRight(lsEnv[i],ψ[i],H[i])
            totaltruncerror = max(totaltruncerror,temptruncerr)
        end
        println(">>>>>> finished >>>>>>")

        println("<<<<<< begin <<<<<<")
        for i in L:-1:2
            Eg,Ev = LocalEigen(H[i-1:i],lsEnv[i-1],lsEnv[i+1],LanczosLevel)
            ψ[i-1:i],temptruncerr = collect(LeftSVD(Ev,D_MPS))
            #Eg, ψ[i-1:i] = LeftUpdateDMRG1(ψ[i-1],H[i],lsEnv[i],lsEnv[i+1],LanczosLevel,D_MPS)
            lsEnv[i] = PushLeft(lsEnv[i+1],ψ[i],H[i])
            totaltruncerror = max(totaltruncerror,temptruncerr)
        end
        println("<<<<<< finished <<<<<<")

        println("sweep $i finished, Eg = $Eg, time consumed $(round(time()-start_time;digits=2)), max truncation error = $(totaltruncerror)")
        lsE[i] = Eg
    end

    return ψ,lsE
end

function LocalEigen(Hi::AbstractTensorMap,
    EnvL::AbstractTensorMap,EnvR::AbstractTensorMap,
    LanczosLevel::Int64)

    effH = EffHam(Hi,EnvL,EnvR)
    return groundEig(effH,LanczosLevel)
end

function LocalEigen(Hi::Vector{AbstractTensorMap},
    EnvL::AbstractTensorMap,EnvR::AbstractTensorMap,
    LanczosLevel::Int64)

    effH = EffHam(Hi,EnvL,EnvR)
    return groundEig(effH,LanczosLevel)
end


function RightUpdateDMRG1(nextψ::AbstractTensorMap,Hi::AbstractTensorMap,
    EnvL::AbstractTensorMap,EnvR::AbstractTensorMap,
    LanczosLevel::Int64,D_MPS::Int64)

    effH = EffHam(Hi,EnvL,EnvR)
    Eg,Ev = groundEig(effH,LanczosLevel)
    MPSs = RightMove(nextψ,Ev,D_MPS)

    return Eg, MPSs
end



function LeftUpdateDMRG1(nextψ::AbstractTensorMap,Hi::AbstractTensorMap,
    EnvL::AbstractTensorMap,EnvR::AbstractTensorMap,
    LanczosLevel::Int64,D_MPS::Int64)

    effH = EffHam(Hi,EnvL,EnvR)
    Eg,Ev = groundEig(effH,LanczosLevel)
    MPSs = LeftMove(nextψ,Ev,D_MPS)

    return Eg, MPSs
end


