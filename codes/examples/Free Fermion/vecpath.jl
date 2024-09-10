using TensorKit,JLD2,LinearAlgebra,FiniteLattices,CairoMakie

include("model.jl")
include("../../src/iMPS.jl")

# apply MPO to MPS
function Apply(MPS::Vector,MPO::Vector)
    

    return finalMPS
end


Latt = YCSqua(3,3)

CKdaggMPO(Latt,[0.2,0.6])


