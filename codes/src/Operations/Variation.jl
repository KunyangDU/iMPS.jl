function VariContract(Opr::Vector{AbstractTensorMap{ComplexSpace,2,2}},ψ₀::Vector,D_MPS::Int64;
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


function VariPlusMPO(
    Opr1::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}},
    Opr2::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}},
    d::Int64,D_MPO::Int64;
    MaxIter::Int64=2
    )

    L = if length(Opr1)==length(Opr2) length(Opr1) else @error "wrong" end

    testOpr = RandMPO(L,d)

    tempOpr1 = deepcopy(Opr1)
    tempOpr2 = deepcopy(Opr2)
    lsEnv1 = vcat(LeftLsEnv(tempOpr1,testOpr,1),RightLsEnv(tempOpr1,testOpr,1))
    lsEnv2 = vcat(LeftLsEnv(tempOpr2,testOpr,1),RightLsEnv(tempOpr2,testOpr,1))

    totaltruncerror = 0
    for iter in 1:MaxIter
        println("sweep $iter")
        start_time = time()

        println(">>>>>> begin >>>>>>")
        for i in 1:L-1
            localOpr = LocalMerge(lsEnv1[i],tempOpr1[i:i+1]...,lsEnv1[i+2]) + LocalMerge(lsEnv2[i],tempOpr2[i:i+1]...,lsEnv2[i+2])
            testOpr[i:i+1],temptruncerr_test = mySVD(localOpr,"right",D_MPO)

            tempOpr1[i:i+1] = Move(tempOpr1[i:i+1]...)
            tempOpr2[i:i+1] = Move(tempOpr2[i:i+1]...)

            lsEnv1[i+1] = PushRight(lsEnv1[i],tempOpr1[i],testOpr[i])
            lsEnv2[i+1] = PushRight(lsEnv2[i],tempOpr2[i],testOpr[i])

            totaltruncerror = max(totaltruncerror,temptruncerr_test)
        end 
        println(">>>>>> finished >>>>>>")        
        
        println("<<<<<< begin <<<<<<")
        for i in L:-1:2
            localOpr = LocalMerge(lsEnv1[i-1],tempOpr1[i-1:i]...,lsEnv1[i+1]) + LocalMerge(lsEnv2[i-1],tempOpr2[i-1:i]...,lsEnv2[i+1])
            testOpr[i-1:i],temptruncerr_test = mySVD(localOpr,"left",D_MPO)

            tempOpr1[i-1:i] = Move(tempOpr1[i-1:i]...)
            tempOpr2[i-1:i] = Move(tempOpr2[i-1:i]...)

            lsEnv1[i] = PushLeft(lsEnv1[i+1],tempOpr1[i],testOpr[i])
            lsEnv2[i] = PushLeft(lsEnv2[i+1],tempOpr2[i],testOpr[i])

            totaltruncerror = max(totaltruncerror,temptruncerr_test)
        end
        println("<<<<<< finished <<<<<<")

        println("sweep $iter finished, time consumed $(round(time()-start_time;digits=2)), max truncation error = $(totaltruncerror)")
    end

    return testOpr
end

function VariProdMPO(
    Opr1::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}},
    Opr2::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}},
    d::Int64,D_MPO::Int64;
    MaxIter::Int64=2
    )

    L = if length(Opr1)==length(Opr2) length(Opr1) else @error "wrong" end

    testOpr = RandMPO(L,d)

    tempOpr1 = deepcopy(Opr1)
    tempOpr2 = deepcopy(Opr2)
    lsEnv = vcat(LeftLsEnv(tempOpr1,tempOpr2,testOpr,1),RightLsEnv(tempOpr1,tempOpr2,testOpr,1))

    totaltruncerror = 0
    for iter in 1:MaxIter
        println("sweep $iter")
        start_time = time()

        println(">>>>>> begin >>>>>>")
        for i in 1:L-1
            localOpr = LocalMerge(lsEnv[i],tempOpr1[i:i+1]...,tempOpr2[i:i+1]...,lsEnv[i+2])
            testOpr[i:i+1],temptruncerr_test = mySVD(localOpr,"right",D_MPO)

            tempOpr1[i:i+1] = Move(tempOpr1[i:i+1]...)
            tempOpr2[i:i+1] = Move(tempOpr2[i:i+1]...)

            lsEnv[i+1] = PushRight(lsEnv[i],tempOpr1[i],tempOpr2[i],testOpr[i])

            totaltruncerror = max(totaltruncerror,temptruncerr_test)
        end 
        println(">>>>>> finished >>>>>>")        
        
        println("<<<<<< begin <<<<<<")
        for i in L:-1:2
            localOpr = LocalMerge(lsEnv[i-1],tempOpr1[i-1:i]...,tempOpr2[i-1:i]...,lsEnv[i+1])
            testOpr[i-1:i],temptruncerr_test = mySVD(localOpr,"left",D_MPO)

            tempOpr1[i-1:i] = Move(tempOpr1[i-1:i]...)
            tempOpr2[i-1:i] = Move(tempOpr2[i-1:i]...)

            lsEnv[i] = PushLeft(lsEnv[i+1],tempOpr1[i],tempOpr2[i],testOpr[i])

            totaltruncerror = max(totaltruncerror,temptruncerr_test)
        end
        println("<<<<<< finished <<<<<<")

        println("sweep $iter finished, time consumed $(round(time()-start_time;digits=2)), max truncation error = $(totaltruncerror)")
    end

    return testOpr
    
end

function VariPlusMPS(
    ψ1::Vector{Union{AbstractTensorMap{ComplexSpace,0,3},AbstractTensorMap{ComplexSpace,1,2}}},
    ψ2::Vector{Union{AbstractTensorMap{ComplexSpace,0,3},AbstractTensorMap{ComplexSpace,1,2}}},
    d::Int64,D_MPS::Int64;
    MaxIter::Int64 = 2
    )

    L = if length(ψ1)==length(ψ2) length(ψ1) else @error "wrong" end

    testψ = RandMPS(L;d=d)

    tempψ1 = deepcopy(ψ1)
    tempψ2 = deepcopy(ψ2)
    lsEnv1 = vcat(LeftLsEnv(tempψ1,testψ,1),RightLsEnv(tempψ1,testψ,1))
    lsEnv2 = vcat(LeftLsEnv(tempψ2,testψ,1),RightLsEnv(tempψ2,testψ,1))

    totaltruncerror = 0
    for iter in 1:MaxIter
        println("sweep $iter")
        start_time = time()

        println(">>>>>> begin >>>>>>")
        for i in 1:L-1
            localψ = LocalMerge(lsEnv1[i],tempψ1[i:i+1]...,lsEnv1[i+2]) + LocalMerge(lsEnv2[i],tempψ2[i:i+1]...,lsEnv2[i+2])
            testψ[i:i+1],temptruncerr_test = mySVD(localψ,"right",D_MPS)

            tempψ1[i:i+1] = Move(tempψ1[i:i+1]...)
            tempψ2[i:i+1] = Move(tempψ2[i:i+1]...)

            lsEnv1[i+1] = PushRight(lsEnv1[i],tempψ1[i],testψ[i])
            lsEnv2[i+1] = PushRight(lsEnv2[i],tempψ2[i],testψ[i])

            totaltruncerror = max(totaltruncerror,temptruncerr_test)
        end 
        println(">>>>>> finished >>>>>>")
        
        println("<<<<<< begin <<<<<<")
        for i in L:-1:2
            localOpr = LocalMerge(lsEnv1[i-1],tempψ1[i-1:i]...,lsEnv1[i+1]) + LocalMerge(lsEnv2[i-1],tempψ2[i-1:i]...,lsEnv2[i+1])
            testψ[i-1:i],temptruncerr_test = mySVD(localOpr,"left",D_MPO)

            tempψ1[i-1:i] = Move(tempψ1[i-1:i]...)
            tempψ2[i-1:i] = Move(tempψ2[i-1:i]...)

            lsEnv1[i] = PushLeft(lsEnv1[i+1],tempψ1[i],testψ[i])
            lsEnv2[i] = PushLeft(lsEnv2[i+1],tempψ2[i],testψ[i])

            totaltruncerror = max(totaltruncerror,temptruncerr_test)
        end
        println("<<<<<< finished <<<<<<")

        println("sweep $iter finished, time consumed $(round(time()-start_time;digits=2)), max truncation error = $(totaltruncerror)")
    end
    return testψ
end
