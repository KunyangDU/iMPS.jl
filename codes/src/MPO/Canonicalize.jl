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