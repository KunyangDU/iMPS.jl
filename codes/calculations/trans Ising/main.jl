using TensorKit,LinearAlgebra,JLD2

include("../iMPS/iMPS.jl")
include("model.jl")

L = 12

d = 2
D_MPO = 3
D_MPS = 2^4

J = -1.0

for h in -2.0:0.2:2.0
    
    H = IsingHam(L;J=J,h=h)
    ψ = initialMPS(L,d,D_MPS)

    LanczosLevel = 16
    Nsweep = 5

    ψ,lsE = sweepDMRG1(ψ,H,Nsweep,LanczosLevel,D_MPS)

    @save "trans Ising/data/tψ_D=$(D_MPS)_L=$(L)_J=$(J)_h=$(h).jld2" ψ
    @save "trans Ising/data/tlsE_D=$(D_MPS)_L=$(L)_J=$(J)_h=$(h).jld2" lsE
end


