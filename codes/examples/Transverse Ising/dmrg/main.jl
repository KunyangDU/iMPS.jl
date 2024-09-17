using JLD2,TensorKit
include("../../../src/iMPS.jl")
include("../model.jl")

L = 12

d = 2
D_MPO = 3
D_MPS = 2^3

J = -1.0

for h in 0.0
    
    H = HamMPO(L;J=J,h=h)
    ψ = RandMPS(L)

    LanczosLevel = 16
    Nsweep = 3

    ψ,lsE = sweepDMRG2(ψ,H,Nsweep,LanczosLevel,D_MPS)

    showQuantSweep(lsE;name="Eg sweep")

    #@save "examples/Transverse Ising/data/dmrg/ψ_D=$(D_MPS)_L=$(L)_J=$(J)_h=$(h).jld2" ψ
    #@save "examples/Transverse Ising/data/dmrg/lsE_D=$(D_MPS)_L=$(L)_J=$(J)_h=$(h).jld2" lsE
end


