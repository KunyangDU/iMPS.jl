using TensorKit,LinearAlgebra,JLD2

include("../iMPS/iMPS.jl")
include("model.jl")

L = 12
params = (J = 0.0,h = 0.5)
D_MPS = 128

d = 2
D_magmom = 2

ψ = load("trans Ising/data/ψ_D=$(D_MPS)_L=$(L)_J=$(params[1])_h=$(params[2]).jld2")["ψ"]
MM = IsingMagmom(L,d,D_magmom)
reduM = reduHam(ψ,MM,1)

@show reduM' ≈ reduM

M = @tensor ψ[1][2,1]*reduM[1,2,3,4]*ψ[1]'[4,3]


