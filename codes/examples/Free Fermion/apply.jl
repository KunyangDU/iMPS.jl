using TensorKit,JLD2,LinearAlgebra,FiniteLattices,CairoMakie

include("model.jl")
include("../../src/iMPS.jl")

# apply MPO to MPS
#= function VariContract(MPO::Vector,MPS::Vector,D_MPS::Int64;
    d::Number=2,MaxIter::Int64=5)

    L = length(MPS)
    ψtest = RandMPS(L;d=d)

    # calculate the Environment
    lsEnv = vcat(LeftLsEnv(MPS,MPO,ψtest,1),RightLsEnv(MPS,MPO,ψtest,1))

    for iter in 1:MaxIter
        println("sweep $iter")
        start_time = time()

        println(">>>>>> begin >>>>>>")
        for i in 1:L-1
            Ev = Apply(LocalMerge(MPS[i:i+1]...),EffHam(MPO[i:i+1],lsEnv[i],lsEnv[i+2]))
            ψtest[i:i+1],temptruncerr_test = RightSVD(Ev,D_MPS)
            lsEnv[i+1] = PushRight(lsEnv[i],MPS[i],MPO[i],ψtest[i])
            MPS[i:i+1],temptruncerr_MPS = RightMove(MPS[i+1],Ev,D_MPS)
            totaltruncerror = max(totaltruncerror,temptruncerr_test,temptruncerr_MPS)
        end
        println(">>>>>> finished >>>>>>")

        
        println("<<<<<< begin <<<<<<")
        for i in L:-1:2
            Ev = Apply(LocalMerge(MPS[i-1:i]...),EffHam(MPO[i-1:i],lsEnv[i-1],lsEnv[i+1]))
            ψ[i-1:i],temptruncerr_test = LeftSVD(Ev,D_MPS)
            lsEnv[i] = PushLeft(lsEnv[i+1],MPS[i],MPO[i],ψtest[i])
            MPS[i-1:i],temptruncerr_MPS = LeftMove(MPS[i-1],Ev,D_MPS)
            totaltruncerror = max(totaltruncerror,temptruncerr_test,temptruncerr_MPS)
        end
        println("<<<<<< finished <<<<<<")

        println("sweep $iter finished, time consumed $(round(time()-start_time;digits=2)), max truncation error = $(totaltruncerror)")
    end

    return ψtest
end =#








