function canonicalize(Opr::Vector{AbstractTensorMap{ComplexSpace,2,2}})
    #= 
    canonicalize!
    从SETTN开始所有的MPO都是正则化的，包括RandMPO
    其它零温计算的MPO都不正则化
    =#
    Opr = convert(Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}},Opr)
    Opr[end] = Move(Opr[end],InitialEnv(Vector{AbstractTensorMap{ComplexSpace,1,3}}(undef,1)))[1]
    for i in length(Opr):-1:2
        Opr[i-1:i] = Move(Opr[i-1:i]...)
    end

    return Opr
end

function compress(Opr::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}};trunc::Float64 = 1e-3)
    D_MPO = 0
    L = length(Opr)

    for iL in 1:L-1
        Oprm = Contract(Opr[iL:iL+1]...)
        Opr[iL:iL+1],temp_D_MPO = mySVD(Oprm,"right",trunc)
        D_MPO = max(D_MPO,temp_D_MPO)
    end

    for iL in L:-1:2
        Oprm = Contract(Opr[iL-1:iL]...)
        Opr[iL-1:iL],temp_D_MPO = mySVD(Oprm,"left",trunc)
        D_MPO = max(D_MPO,temp_D_MPO)
    end

    return Opr,D_MPO
end