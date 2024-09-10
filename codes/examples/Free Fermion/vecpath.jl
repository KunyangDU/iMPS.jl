using TensorKit,JLD2,LinearAlgebra,FiniteLattices,CairoMakie

include("model.jl")
include("../../src/iMPS.jl")

# apply MPO to MPS
function Apply(MPS::Vector,MPO::Vector)
    

    return finalMPS
end

function test(a)
    return [a,a],a
end


a = 1
b = zeros(3)
c = 0
d = 0
e = 0

c,d,e = test(a)


