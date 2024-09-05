using TensorKit,LinearAlgebra,JLD2

include("../iMPS/iMPS.jl")
include("model.jl")

D_MPS = 128
L = 12

D_magmom = 2
d = 2

lsParams = [(J = 0.0,h = 0.5),(J = 1.0,h = 0.)]
params = lsParams[1]

ψ = load("trans Ising/data/ψ_D=$(D_MPS)_L=$(L)_J=$(params[1])_h=$(params[2]).jld2")["ψ"]
MM = IsingMagmom(L,d,D_magmom)

effM = reduHam(ψ,MM,1)

M = @tensor ψ[1][2,1]*effM[1,2,3,4]*ψ[1]'[4,3]


cψ = []
for mps in ψ
    push!(cψ,mps*mps')
end
norm(cψ[1])
