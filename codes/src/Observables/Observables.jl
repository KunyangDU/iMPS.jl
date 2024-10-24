


function GreenFuncFreq(
    ψ::Vector{Union{AbstractTensorMap{ComplexSpace,1,2},AbstractTensorMap{ComplexSpace,0,3}}},
    Opr::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}},
    Oprdagg::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}},
    H::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}},
    E0::Number,D_MPS::Int64,LanczosLevel::Int64,lsω::Vector;
    τ::Number=1e-2,MaxIter::Int64=200,TruncErr::Number=1e-3,ϵ::Int64 = 1
    )
    lsGF₊ = CorrFuncFreq(ψ,Oprdagg,H,E0,D_MPS,LanczosLevel,lsω; τ=τ,MaxIter=MaxIter,TruncErr=TruncErr,TR=true)
    lsGF₋ = CorrFuncFreq(ψ,Opr,H,E0,D_MPS,LanczosLevel,-lsω; τ=-τ,MaxIter=MaxIter,TruncErr=TruncErr,TR=true)
    lsGF = @. -1im * (lsGF₊ - ϵ * lsGF₋) / 2
    # for the time integral of GF is only half that of correlation function
    # only suitable for spectral function, for the imaginary part is discarded in CorrFuncFreq (make another one!)
    return lsGF
end

function SpecFuncFreq(
    ψ::Vector{Union{AbstractTensorMap{ComplexSpace,1,2},AbstractTensorMap{ComplexSpace,0,3}}},
    Opr::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}},
    Oprdagg::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}},
    H::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}},
    E0::Number,D_MPS::Int64,LanczosLevel::Int64,lsω::Vector;
    τ::Number=1e-2,MaxIter::Int64=200,TruncErr::Number=1e-3,ϵ::Int64 = 1
    )
    lsSF = (-1/pi) .* imag.( GreenFuncFreq(ψ,Opr,Oprdagg,H,E0,D_MPS,LanczosLevel,lsω; τ=-τ,MaxIter=MaxIter,TruncErr=TruncErr,ϵ=ϵ))
    return lsSF
end

function Fidelity(
    ψ::Vector{Union{AbstractTensorMap{ComplexSpace,1,2},AbstractTensorMap{ComplexSpace,0,3}}},
    H::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}},
    E0::Number,D_MPS::Int64,LanczosLevel::Int64;
    τ::Number=1e-2,MaxIter::Int64=200,TruncErr::Number=1e-3,
    )

    L = length(H)

    lsF = Vector{ComplexF64}(undef,1)
    lsF[1] = InnerProd(ψ,ψ) / 2

    lst = Vector{Float64}(undef,1)
    lst[1] = 0.0
    
    ψ₀ = deepcopy(ψ)

    lsEnv = vcat(LeftLsEnv(ψ,H,1),RightLsEnv(ψ,H,1))

    totaltruncerror = 0
    for iter in 1:MaxIter
        start_time = time()

        println("evolution $iter")
        totaltruncerror = TDVP2!(ψ,H,lsEnv,τ,D_MPS,LanczosLevel,totaltruncerror)        
        println("evolution $iter finished, time consumed $(round(time()-start_time;digits=2))s, max truncation error = $(totaltruncerror)")

        totaltruncerror > TruncErr && break

        push!(lst,lst[end] + τ)
        push!(lsF,InnerProd(ψ₀,ψ)*exp(-1im*E0*lst[end]))
    end
    return lsF,lst
end

function CorrFuncTime(
    ψ::Vector{Union{AbstractTensorMap{ComplexSpace,1,2},AbstractTensorMap{ComplexSpace,0,3}}},
    Opr::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}},
    H::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}},
    E0::Number,D_MPS::Int64,LanczosLevel::Int64;
    τ::Number=1e-2,MaxIter::Int64=200,TruncErr::Number=1e-3,
    )
    Oprψ = VariContract(Opr,ψ,D_MPS)
    return Fidelity(Oprψ,H,E0,D_MPS,LanczosLevel; τ=τ,MaxIter=MaxIter,TruncErr=TruncErr)
end

function CorrFuncFreq(
    ψ::Vector{Union{AbstractTensorMap{ComplexSpace,1,2},AbstractTensorMap{ComplexSpace,0,3}}},
    Opr::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}},
    H::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}},
    E0::Number,D_MPS::Int64,LanczosLevel::Int64,lsω::Vector;
    τ::Number=1e-2,MaxIter::Int64=200,TruncErr::Number=1e-3,TR::Bool = false
    )
    if !TR
        lsCF₊ = Time2Freq(CorrFuncTime(ψ,Opr,H,E0,D_MPS,LanczosLevel; τ=τ,MaxIter=MaxIter,TruncErr=TruncErr)...,lsω)
        lsCF₋ = Time2Freq(CorrFuncTime(ψ,Opr,H,E0,D_MPS,LanczosLevel; τ=-τ,MaxIter=MaxIter,TruncErr=TruncErr)...,lsω)
        lsCF = lsCF₊ + lsCF₋
    else
        lsCF = 2 .* real.(Time2Freq(CorrFuncTime(ψ,Opr,H,E0,D_MPS,LanczosLevel; τ=τ,MaxIter=MaxIter,TruncErr=TruncErr)...,lsω))
    end

    return lsCF
end

function Time2Freq(lsF::Vector,lst::Vector,lsω::Vector)
    τ = lst[2] - lst[1]
    return [sum(@.(lsF.*lst*exp(-1im*ω*lst)))*τ for ω in lsω]
end



#= 
function GreenFuncTDVP2(ψ::Vector,H::Vector,τ::Number,
    TruncErr::Number,MaxIter::Int64,D_MPS::Int64,LanczosLevel::Int64)
    # why truncerr is always 0?
    # calculate the ⟨ψ| exp(-iHt) |ψ⟩

    L = length(H)
    
    lsGt = Vector{ComplexF64}(undef,1)
    lsGt[1] = InnerProd(ψ,ψ) / 2
    lst = Vector{Float64}(undef,1)
    lst[1] = 0.0

    
    ψ₀ = deepcopy(ψ)

    lsEnv = vcat(LeftLsEnv(ψ,H,1),RightLsEnv(ψ,H,1))

    
    totaltruncerror = 0
    for iter in 1:MaxIter

        start_time = time()
        println("evolution $iter")

        println(">>>>>> begin >>>>>>")
        for i in 1:L-1
            if iter != 1 && i==1
                ψ[i:i+1],temptruncerr = RightUpdateTDVP2(ψ[i:i+1],H[i:i+1],lsEnv[i],lsEnv[i+2],2*τ/2,D_MPS,LanczosLevel;τback=τ/2)
            else
                ψ[i:i+1],temptruncerr = RightUpdateTDVP2(ψ[i:i+1],H[i:i+1],lsEnv[i],lsEnv[i+2],τ/2,D_MPS,LanczosLevel)
            end
            lsEnv[i+1] = PushRight(lsEnv[i],ψ[i],H[i])
            totaltruncerror = max(totaltruncerror,temptruncerr)
        end
        println(">>>>>> finished >>>>>>")

        println("<<<<<< begin <<<<<<")
        for i in L:-1:2
            if i == L
                ψ[i-1:i],temptruncerr = LeftUpdateTDVP2(ψ[i-1:i],H[i-1:i],lsEnv[i-1],lsEnv[i+1],2*τ/2,D_MPS,LanczosLevel;τback=τ/2)
            else
                ψ[i-1:i],temptruncerr = LeftUpdateTDVP2(ψ[i-1:i],H[i-1:i],lsEnv[i-1],lsEnv[i+1],τ/2,D_MPS,LanczosLevel)
            end
            lsEnv[i] = PushLeft(lsEnv[i+1],ψ[i],H[i])
            totaltruncerror = max(totaltruncerror,temptruncerr)
        end
        println("<<<<<< finished <<<<<<")

        println("evolution $iter finished, time consumed $(round(time()-start_time;digits=2))s, max truncation error = $(totaltruncerror)")
        totaltruncerror > TruncErr && break
        push!(lsGt,InnerProd(ψ₀,ψ))
        push!(lst,lst[end] + τ)
        # without a -iexp(...)
    end

    return lsGt,lst
end
=#

