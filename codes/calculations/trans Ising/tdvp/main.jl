using TensorKit,LinearAlgebra,JLD2

include("../../iMPS/iMPS.jl")
include("../model.jl")

L = 12

d = 2
D_MPO = 3
D_MPS = 2^4

J = -1.0


