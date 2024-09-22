

function MonoTrace(Opr1::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}})
    
    EnvR = RightEnv(Opr1,1)
    EnvL = LeftEnv(Opr1,1)

    tr = @tensor EnvL[1]*Opr1[1][2,1,2,3]*EnvR[3]

    return ApproxReal(tr[1])
end

function Trace(Opr1::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}})
    
    EnvR = RightEnv(Opr1,Opr1,1)
    EnvL = LeftEnv(Opr1,Opr1,1)

    tr = @tensor EnvL[1,2]*Opr1[1][4,1,3,5]*Opr1[1]'[2,3,6,4]*EnvR[5,6]

    return ApproxReal(tr[1] / norm(tr[1])) * norm(tr[1])
end

function Trace(
    Opr1::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}},
    Opr2::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}};
    Z::Number=1
    )

    EnvR = RightEnv(Opr1,Opr2,Opr1,1)
    EnvL = LeftEnv(Opr1,Opr2,Opr1,1)

    tr = @tensor EnvL[1,2,4]*Opr1[1][6,1,3,7]*Opr2[1][3,2,5,8]*Opr1[1]'[4,5,9,6]*EnvR[7,8,9]

    return ApproxReal(tr[1]/Z)
end


