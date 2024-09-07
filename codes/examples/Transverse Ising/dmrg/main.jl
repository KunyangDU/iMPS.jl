using TensorKit,LinearAlgebra,JLD2

include("../../../src/iMPS.jl")
include("../model.jl")

L = 12

d = 2
D_MPO = 3
D_MPS = 2^5

J = 1.0

for h in -2.0:0.1:2.0
    
    H = HamMPO(L;J=J,h=h)
    ψ = RandMPS(L,d)

    LanczosLevel = 16
    Nsweep = 3

    ψ,lsE1 = sweepDMRG2(ψ,H,Nsweep,LanczosLevel,D_MPS)
    ψ,lsE2 = sweepDMRG1(ψ,H,Nsweep,LanczosLevel,D_MPS)

    showQuantSweep(lsE2;name="Eg sweep2")
    showQuantSweep(lsE1;name="Eg sweep1")
    lsE = vcat(lsE1,lsE2)

    @save "examples/Transverse Ising/data/dmrg/ψ_D=$(D_MPS)_L=$(L)_J=$(J)_h=$(h).jld2" ψ
    @save "examples/Transverse Ising/data/dmrg/lsE_D=$(D_MPS)_L=$(L)_J=$(J)_h=$(h).jld2" lsE
end


