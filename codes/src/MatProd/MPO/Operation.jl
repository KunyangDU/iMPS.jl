
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


function WeightProd(MPS::Vector,weight::Number)
    return (vcat([weight],ones(length(MPS)-1)) .* MPS)
end

function WeightSum(MPSs::Vector,weights::Vector,D_MPS::Int64)
    d = dims(domain(MPSs[1][1]))


    MPS = VariPlusMPS(WeightProd(MPSs[1],weights[1]), WeightProd(MPSs[2],weights[2]),d,D_MPS)

    for i in eachindex(MPSs)[3:end]
        MPS = VariPlusMPS(MPS,WeightProd(MPSs[i],weights[i]),d,D_MPS)
    end
    
    return MPS
end


