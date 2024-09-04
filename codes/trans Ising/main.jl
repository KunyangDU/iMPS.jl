using TensorKit,LinearAlgebra

include("../iMPS/iMPS.jl")
include("model.jl")

L = 12
J = 1.0
h = 0.5
d = 2
D_MPO = 3
D_MPS = 128

H = tranIsingH(L,J,h,d,D_MPO)
ψ = initialMPS(L,d,D_MPS)

LanczosLevel = 8
Nsweep = 10

ψ,lsE = sweep1(ψ,H,Nsweep,LanczosLevel,D_MPS)

for i in eachindex(lsE)
    println(lsE[i])
end

@save "trans Ising/data/ψ_D=$(D_MPS)_L=$(L).jld2" ψ


