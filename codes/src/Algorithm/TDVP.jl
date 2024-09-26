
function MixTDVP(ψ::Vector,H::Vector,
    t::Number,Nt::Int64,
    D_MPS::Int64;prop::NTuple{2,Int64} = (1,5),TruncErr::Number = 1e-5)

    N2, N1 = convert.(Int64,Nt .* prop ./ sum(prop))
    t2, t1 = convert.(Float64,t .* prop ./ sum(prop))

    lsψ2,lst2 = sweepTDVP2(ψ,H,t2,N2,D_MPS;TruncErr=TruncErr)
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
    totaltruncerror = 0
    for iNt in 2:Nt

        start_time = time()
        println("evolution $iNt, t = $(round(lst[iNt];digits=3))/J")

        println(">>>>>> begin >>>>>>")
        for i in 1:L-1
            if iNt != 1 && i == 1
                ψ[i:i+1],temptruncerr = RightUpdateTDVP1(ψ[i:i+1],H[i],lsEnv[i],lsEnv[i+1],2*τ,D_MPS;τback = τ)
            else
                ψ[i:i+1],temptruncerr = RightUpdateTDVP1(ψ[i:i+1],H[i],lsEnv[i],lsEnv[i+1],τ,D_MPS)
            end
            lsEnv[i+1] = PushRight(lsEnv[i],ψ[i],H[i])

            totaltruncerror = max(totaltruncerror,temptruncerr)
        end
        println(">>>>>> finished >>>>>>")

        println("<<<<<< begin <<<<<<")
        for i in L:-1:2
            if i == L 
                ψ[i-1:i],temptruncerr = LeftUpdateTDVP1(ψ[i-1:i],H[i],lsEnv[i],lsEnv[i+1],2*τ,D_MPS;τback = τ)
            else
                ψ[i-1:i],temptruncerr = LeftUpdateTDVP1(ψ[i-1:i],H[i],lsEnv[i],lsEnv[i+1],τ,D_MPS)
            end
            lsEnv[i] = PushLeft(lsEnv[i+1],ψ[i],H[i])

            totaltruncerror = max(totaltruncerror,temptruncerr)
        end
        println("<<<<<< finished <<<<<<")

        println("evolution $iNt finished, time consumed $(round(time()-start_time;digits=2))s, max truncation error = $(totaltruncerror)")

        lsψ[iNt] = deepcopy(ψ)
    end

    return lsψ,lst
end

function RightUpdateTDVP1(ψs::Vector,Hi::AbstractTensorMap,
    EnvL::AbstractTensorMap,EnvR::AbstractTensorMap,τ::Number,
    D_MPS::Int64;τback::Number = τ)

    effH = EffHam(Hi,EnvL,EnvR)
    Aτ = Apply(ψs[1],EvolveOpr(effH,τ))

    MPSs,truncerr = RightSVD(Aτ,D_MPS)
    thisMPS,Στ = MPSs

    effH0 = EffHam(PushRight(EnvL,thisMPS,Hi),EnvR)
    Σ = Apply(Στ,EvolveOpr(effH0,-τback))

    return RightMerge(Σ,thisMPS,ψs[2]),truncerr
end


function LeftUpdateTDVP1(ψs::Vector,Hi::AbstractTensorMap,
    EnvL::AbstractTensorMap,EnvR::AbstractTensorMap,τ::Number,
    D_MPS::Int64;τback::Number = τ)

    effH = EffHam(Hi,EnvL,EnvR)
    Aτ = Apply(ψs[2],EvolveOpr(effH,τ))

    MPSs,truncerr = LeftSVD(Aτ,D_MPS)
    Στ,thisMPS = MPSs

    effH0 = EffHam(EnvL,PushLeft(EnvR,thisMPS,Hi))
    Σ = Apply(Στ,EvolveOpr(effH0,-τback))

    return LeftMerge(Σ,thisMPS,ψs[1]),truncerr
end
############ NORMALIZED ##################

function TDVP2!(
    ψ::Vector{Union{AbstractTensorMap{ComplexSpace,1,2},AbstractTensorMap{ComplexSpace,0,3}}},
    H::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}},
    t::Number,Nt::Int64,
    D_MPS::Int64,LanczosLevel::Int64;TruncErr::Number=1e-3)

    L = length(H)
    
    lsψ = Vector{Vector}(undef,1)
    lst = Vector{Float64}(undef,1)
    τ = t/(Nt-1)

    lsψ[1] = deepcopy(ψ)
    lst[1] = 0.0

    lsEnv = vcat(LeftLsEnv(ψ,H,1),RightLsEnv(ψ,H,1))

    totaltruncerror = 0
    for iNt in 2:Nt

        start_time = time()
        println("evolution $iNt, t = $(round(lst[iNt-1]+τ;digits=3))/J")
        totaltruncerror = TDVP2!(ψ,H,lsEnv,τ,D_MPS,LanczosLevel,totaltruncerror)        
        println("evolution $iNt finished, time consumed $(round(time()-start_time;digits=2))s, max truncation error = $(totaltruncerror)")

        totaltruncerror > TruncErr && break
        push!(lsψ,deepcopy(ψ))
        push!(lst,lst[end] + τ)
    end

    return lsψ,lst
end


function TDVP2!(
    ψ::Vector{Union{AbstractTensorMap{ComplexSpace,1,2},AbstractTensorMap{ComplexSpace,0,3}}},
    H::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}},
    lsEnv::Vector{AbstractTensorMap{ComplexSpace,2,1}},
    τ::Number,D_MPS::Int64,LanczosLevel::Int64,totaltruncerror::Number)
    
    L = length(ψ)
    println(">>>>>> begin >>>>>>")
    for i in 1:L-1
        ψ[i:i+1],H[i:i+1],temptruncerr = RightUpdateTDVP2(ψ[i:i+1],H[i:i+1],lsEnv[i],lsEnv[i+2],τ / 2,D_MPS,LanczosLevel)
        lsEnv[i+1] = PushRight(lsEnv[i],ψ[i],H[i])
        totaltruncerror = max(totaltruncerror,temptruncerr)
    end
    ψ[L] = Evolve(ψ[L],H[L],lsEnv[L],lsEnv[L+1],τ / 2,LanczosLevel)
    println(">>>>>> finished >>>>>>")

    println("<<<<<< begin <<<<<<")
    for i in L:-1:2
        ψ[i-1:i],H[i-1:i],temptruncerr = LeftUpdateTDVP2(ψ[i-1:i],H[i-1:i],lsEnv[i-1],lsEnv[i+1],τ / 2,D_MPS,LanczosLevel)
        lsEnv[i] = PushLeft(lsEnv[i+1],ψ[i],H[i])
        totaltruncerror = max(totaltruncerror,temptruncerr)
    end
    ψ[1] = Evolve(ψ[1],H[1],lsEnv[1],lsEnv[2],τ / 2,LanczosLevel)
    println("<<<<<< finished <<<<<<")

    return totaltruncerror 
end


function RightUpdateTDVP2(ψs::Vector,Hi::Vector,
    EnvL::AbstractTensorMap,EnvR::AbstractTensorMap,τ::Number,
    D_MPS::Int64,LanczosLevel::Int64;
    τback::Number = τ)

    Aτ = Evolve(Contract(ψs...),Hi,EnvL,EnvR,τ,LanczosLevel)
    MPSs,truncerr = mySVD(Aτ,"right",D_MPS)
    thisMPS,Στ = MPSs 
    Hi = Move(Hi...)
    Σ = Evolve(Στ,Hi[2],PushRight(EnvL,thisMPS,Hi[1]),EnvR,-τback,LanczosLevel)

    return [thisMPS,Σ],Hi,truncerr
end

function LeftUpdateTDVP2(ψs::Vector,Hi::Vector,
    EnvL::AbstractTensorMap,EnvR::AbstractTensorMap,τ::Number,
    D_MPS::Int64,LanczosLevel::Int64;
    τback::Number = τ)

    Aτ = Evolve(Contract(ψs...),Hi,EnvL,EnvR,τ,LanczosLevel)
    MPSs,truncerr = LeftSVD(Aτ,D_MPS)
    Στ,thisMPS = MPSs
    Hi = Move(Hi...)
    Σ = Evolve(Στ,Hi[1],EnvL,PushLeft(EnvR,thisMPS,Hi[2]),-τback,LanczosLevel)

    return [Σ,thisMPS],Hi,truncerr
end


