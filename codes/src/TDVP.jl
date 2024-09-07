
function MixTDVP(ψ::Vector,H::Vector,
    t::Number,Nt::Int64,
    D_MPS::Int64;prop::NTuple{2,Int64} = (1,5))

    N2, N1 = convert.(Int64,Nt .* prop ./ sum(prop))
    t2, t1 = convert.(Float64,t .* prop ./ sum(prop))

    lsψ2,lst2 = sweepTDVP2(ψ,H,t2,N2,D_MPS)
    lsψ1,lst1 = sweepTDVP1(deepcopy(lsψ2[end]),H,t1,N1,D_MPS)
    lsψ = vcat(lsψ2,lsψ1)
    lst = vcat(lst2,lst2[end] .+ lst1)
    return lsψ,lst
end

function sweepTDVP1(ψ::Vector,H::Vector,
    t::Number,Nt::Int64,
    D_MPS::Int64)

    L = length(H)
    
    lsψ = Vector{Vector}(undef,Nt)
    lst = collect(range(0,t,Nt))
    τ = t/(Nt-1)/2

    lsψ[1] = deepcopy(ψ)

    lsEnv = vcat(LeftLsEnv(ψ,H,1),RightLsEnv(ψ,H,1))

    for iNt in 2:Nt
        start_time = time()
        println("evolution $iNt, t = $(round(lst[iNt];digits=3))/J")

        println(">>>>>> begin >>>>>>")
        for i in 1:L-1
            if iNt != 1 && i == 1
                ψ[i:i+1] .= RightUpdateTDVP1(ψ[i:i+1],H[i],lsEnv[i],lsEnv[i+1],2*τ,D_MPS;τback = τ)
            else
                ψ[i:i+1] .= RightUpdateTDVP1(ψ[i:i+1],H[i],lsEnv[i],lsEnv[i+1],τ,D_MPS)
            end
            lsEnv[i+1] = PushRight(lsEnv[i],ψ[i],H[i])
        end
        println(">>>>>> finished >>>>>>")

        println("<<<<<< begin <<<<<<")
        for i in L:-1:2
            if i == L 
                ψ[i-1:i] .= LeftUpdateTDVP1(ψ[i-1:i],H[i],lsEnv[i],lsEnv[i+1],2*τ,D_MPS;τback = τ)
            else
                ψ[i-1:i] .= LeftUpdateTDVP1(ψ[i-1:i],H[i],lsEnv[i],lsEnv[i+1],τ,D_MPS)
            end
            lsEnv[i] = PushLeft(lsEnv[i+1],ψ[i],H[i])
        end
        println("<<<<<< finished <<<<<<")

        println("evolution $iNt finished, time consumed $(round(time()-start_time;digits=2))s")

        lsψ[iNt] = deepcopy(ψ)
    end

    return lsψ,lst
end


function sweepTDVP2(ψ::Vector,H::Vector,
    t::Number,Nt::Int64,
    D_MPS::Int64)
    
    lsψ = Vector{Vector}(undef,Nt)
    lst = collect(range(0,t,Nt))
    τ = t/(Nt-1)/2

    lsψ[1] = deepcopy(ψ)

    lsEnv = vcat(LeftLsEnv(ψ,H,1),RightLsEnv(ψ,H,1))

    for iNt in 2:Nt
#=         lsMi = Quant1(deepcopy(ψ), LocalMagmomMPO(),D_MPS)
        @show lsMi =#
        start_time = time()
        println("evolution $iNt, t = $(round(lst[iNt];digits=3))/J")

        println(">>>>>> begin >>>>>>")
        for i in 1:length(H)-1
            if iNt != 1 && i==1
                ψ[i:i+1] .= RightUpdateTDVP2(ψ[i:i+1],H[i:i+1],lsEnv[i],lsEnv[i+2],2*τ,D_MPS;τback=τ)
            else
                ψ[i:i+1] .= RightUpdateTDVP2(ψ[i:i+1],H[i:i+1],lsEnv[i],lsEnv[i+2],τ,D_MPS)
            end
            #ψ[i:i+1] .= RightUpdateTDVP2(ψ[i:i+1],H[i:i+1],lsEnv[i],lsEnv[i+2],τ,D_MPS)
            lsEnv[i+1] = PushRight(lsEnv[i],ψ[i],H[i])
        end
        println(">>>>>> finished >>>>>>")

        println("<<<<<< begin <<<<<<")
        for i in L:-1:2
            if i == L
                ψ[i-1:i] .= LeftUpdateTDVP2(ψ[i-1:i],H[i-1:i],lsEnv[i-1],lsEnv[i+1],2*τ,D_MPS;τback=τ)
            else
                ψ[i-1:i] .= LeftUpdateTDVP2(ψ[i-1:i],H[i-1:i],lsEnv[i-1],lsEnv[i+1],τ,D_MPS)
            end
            #ψ[i-1:i] .= LeftUpdateTDVP2(ψ[i-1:i],H[i-1:i],lsEnv[i-1],lsEnv[i+1],τ,D_MPS)
            lsEnv[i] = PushLeft(lsEnv[i+1],ψ[i],H[i])
        end
#=         ψ[1:2] .= RightUpdateTDVP2(ψ[1:2],H[1:2],lsEnv[1],lsEnv[3],τ,D_MPS)
        lsEnv[2] = PushRight(lsEnv[1],ψ[1],H[1]) =#
        println("<<<<<< finished <<<<<<")

        println("evolution $iNt finished, time consumed $(round(time()-start_time;digits=2))s")

        lsψ[iNt] = deepcopy(ψ)
    end

    return lsψ,lst
end

function RightUpdateTDVP1(ψs::Vector,Hi::AbstractTensorMap,
    EnvL::AbstractTensorMap,EnvR::AbstractTensorMap,τ::Number,
    D_MPS::Int64;τback::Number = τ)

    effH = EffHam(Hi,EnvL,EnvR)
    Aτ = Apply(ψs[1],EvolveOpr(effH,τ))

    thisMPS,Στ = RightSVD(Aτ,D_MPS)

    effH0 = EffHam(PushRight(EnvL,thisMPS,Hi),EnvR)
    Σ = Apply(Στ,EvolveOpr(effH0,-τback))

    return RightMerge(Σ,thisMPS,ψs[2])
end


function LeftUpdateTDVP1(ψs::Vector,Hi::AbstractTensorMap,
    EnvL::AbstractTensorMap,EnvR::AbstractTensorMap,τ::Number,
    D_MPS::Int64;τback::Number = τ)

    effH = EffHam(Hi,EnvL,EnvR)
    Aτ = Apply(ψs[2],EvolveOpr(effH,τ))

    Στ,thisMPS = LeftSVD(Aτ,D_MPS)

    effH0 = EffHam(EnvL,PushLeft(EnvR,thisMPS,Hi))
    Σ = Apply(Στ,EvolveOpr(effH0,-τback))

    return LeftMerge(Σ,thisMPS,ψs[1])
end


function RightUpdateTDVP2(ψs::Vector,Hi::Vector,
    EnvL::AbstractTensorMap,EnvR::AbstractTensorMap,τ::Number,
    D_MPS::Int64;τback::Number = τ)

    ψm = LocalMerge(ψs...)
    effH = EffHam(Hi,EnvL,EnvR)
    Aτ = Apply(ψm,EvolveOpr(effH,τ))

    thisMPS,Στ = RightSVD(Aτ,D_MPS)

    effH1 = EffHam(Hi[2],PushRight(EnvL,thisMPS,Hi[1]),EnvR)
    Σ = Apply(Στ,EvolveOpr(effH1,-τback))

    return thisMPS,Σ
end

function LeftUpdateTDVP2(ψs::Vector,Hi::Vector,
    EnvL::AbstractTensorMap,EnvR::AbstractTensorMap,τ::Number,
    D_MPS::Int64;τback::Number = τ)

    ψm = LocalMerge(ψs...)
    effH = EffHam(Hi,EnvL,EnvR)
    Aτ = Apply(ψm,EvolveOpr(effH,τ))

    Στ,thisMPS = LeftSVD(Aτ,D_MPS)

    effH1 = EffHam(Hi[1],EnvL,PushLeft(EnvR,thisMPS,Hi[2]))
    Σ = Apply(Στ,EvolveOpr(effH1,-τback))

    return Σ,thisMPS
end



#= function LocalEvolve(ψi::AbstractTensorMap,Hi::AbstractTensorMap,
    EnvL::AbstractTensorMap,EnvR::AbstractTensorMap,τ::Number)
    effH = EffHam(Hi,EnvL,EnvR)
    return Apply(ψi,EvolveOpr(effH,τ))
end

function LocalEvolve(ψi::Vector{AbstractTensorMap},
    Hi::Vector{AbstractTensorMap},
    EnvL::AbstractTensorMap,EnvR::AbstractTensorMap,τ::Number)
    
    ψm = LocalMerge(ψi...)
    effH = EffHam(Hi,EnvL,EnvR)

    return Apply(ψm,EvolveOpr(effH,τ))
end

function RightMove(nextψ::AbstractTensor,thisψ::AbstractTensor,
    Hi::AbstractTensorMap,EnvL::AbstractTensorMap,EnvR::AbstractTensorMap,
    τ::Number,D_MPS::Int64)

    Στ,thisMPS = RightSVD(thisψ,D_MPS)

    effH = EffHam(PushRight(EnvL,thisMPS,Hi),EnvR)
    Σ = Apply(Στ,EvolveOpr(effH,-τ))

    return collect(RightMerge(Σ,thisMPS,nextψ))
end

function LeftMove(nextψ::AbstractTensor,thisψ::AbstractTensor,
    Hi::AbstractTensorMap,EnvL::AbstractTensorMap,EnvR::AbstractTensorMap,
    τ::Number,D_MPS::Int64)

    Στ,thisMPS = LeftSVD(thisψ,D_MPS)

    effH = EffHam(EnvL,PushLeft(EnvR,thisMPS,Hi))
    Σ = Apply(Στ,EvolveOpr(effH,-τ))

    return collect(LeftMerge(Σ,thisMPS,nextψ))
end =#

