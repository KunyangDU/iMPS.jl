
function RandMPS(L::Int64;d::Int64=2)
    # suppose center at leftmost
    MPS = Vector{AbstractTensorMap}(undef,L)

    bond = ℂ^1
    phys = (ℂ^d)'

    MPS[1] = Tensor(rand(ComplexF64,1,2,1),bond' ⊗ phys' ⊗ bond') |> x -> permute(x,(),(1,2,3)) / norm(x)
    
    for i in 2:L
        MPS[i] = let
            TensorMap(rand(1,2,1),bond,phys ⊗ bond) |> x -> x / norm(x)
        end
    end
    
    return MPS
end

function ApproxReal(Qi::Number;tol::Float64=1e-5)
    imag(Qi) <= tol && return real(Qi)
    @error "not real"
end


function LocalMerge(ψ1::AbstractTensorMap{ComplexSpace,1,2},
    ψ2::AbstractTensorMap{ComplexSpace,0,3})

    @tensor ψm[-1,-2,-3,-4] ≔ ψ1[1,-1,-2] * ψ2[1,-3,-4]
    return permute(ψm,(),(1,2,3,4))
end

function LocalMerge(ψ1::AbstractTensorMap{ComplexSpace,0,3},
    ψ2::AbstractTensorMap{ComplexSpace,1,2})

    @tensor ψm[-1,-2,-3,-4] ≔ ψ1[-1,-2,1] * ψ2[1,-3,-4]
    return permute(ψm,(),(1,2,3,4))
end


function InnerProd(ψ₁::Vector,ψ₂::Vector)
    EnvL = InitialLeftEnv(;order=2)
    EnvR = RightEnv(ψ₁,ψ₂,1)
    innerprod = @tensor EnvL[1,3]*ψ₁[1][1,5,2]*ψ₂[1]'[3,5,4]*EnvR[2,4]
    return innerprod[1]
end

function VariContract(Opr::Vector,ψ₀::Vector,D_MPS::Int64;
    d::Number=2,MaxIter::Int64=4)

    L = length(Opr)
    ψtest = RandMPS(L;d=d)

    # calculate the Environment
    lsEnv = vcat(LeftLsEnv(ψ₀,Opr,ψtest,1),RightLsEnv(ψ₀,Opr,ψtest,1))
    totaltruncerror = 0
    for iter in 1:MaxIter
        MPS = deepcopy(ψ₀)
        println("sweep $iter")
        start_time = time()

        println(">>>>>> begin >>>>>>")
        for i in 1:L-1
            Ev = Apply(LocalMerge(MPS[i:i+1]...),EffHam(Opr[i:i+1],lsEnv[i],lsEnv[i+2]))
            ψtest[i:i+1],temptruncerr_test = RightSVD(Ev,D_MPS)
            MPS[i:i+1],temptruncerr_MPS = RightMove(MPS[i+1],MPS[i],D_MPS)
            lsEnv[i+1] = PushRight(lsEnv[i],MPS[i],Opr[i],ψtest[i])
            totaltruncerror = max(totaltruncerror,temptruncerr_test,temptruncerr_MPS)
        end 
        println(">>>>>> finished >>>>>>")        
        
        println("<<<<<< begin <<<<<<")
        for i in L:-1:2
            Ev = Apply(LocalMerge(MPS[i-1:i]...),EffHam(Opr[i-1:i],lsEnv[i-1],lsEnv[i+1]))
            ψtest[i-1:i],temptruncerr_test = LeftSVD(Ev,D_MPS)
            MPS[i-1:i],temptruncerr_MPS = LeftMove(MPS[i-1],MPS[i],D_MPS)
            lsEnv[i] = PushLeft(lsEnv[i+1],MPS[i],Opr[i],ψtest[i])
            totaltruncerror = max(totaltruncerror,temptruncerr_test,temptruncerr_MPS)
        end
        println("<<<<<< finished <<<<<<")

        println("sweep $iter finished, time consumed $(round(time()-start_time;digits=2)), max truncation error = $(totaltruncerror)")
    end

    return ψtest
end


function vrange(beginvec::Union{Vector,Tuple},endvec::Union{Vector,Tuple};step::Int64 = 100)
    return hcat([collect(beginvec .+ (endvec .- beginvec) .* t)  for t in range(0,1,step)]...)
end

function vrange(ipath::Matrix;eachstep::Number = 100)
    finalpath = ipath[:,1]

    for ii in 1:size(ipath)[2]-1
        finalpath = hcat(finalpath,vrange(ipath[:,ii],ipath[:,ii+1];step=eachstep+1)[:,2:end])
    end

    return finalpath
end

function pathlength(finalpath::Matrix)
    return cumsum(norm.(eachcol(hcat([0.0,0.0],diff(finalpath,dims = 2)))))
end