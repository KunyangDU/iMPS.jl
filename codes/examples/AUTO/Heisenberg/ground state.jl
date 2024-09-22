include("../../../src/iMPS.jl")
include("model.jl")

Lx = 4
Ly = 1
J = 4
Latt = YCSqua(Lx,Ly)

Nsweep = 5
LanczosLevel = 15
D_MPS = 2^3

H = Hamiltonian(Latt;Jx=J,h=0)

ψ = HeisenbergRandMPS(Latt)

ψ,lsE = sweepDMRG2(ψ,H,Nsweep,LanczosLevel,D_MPS)

showQuantSweep(lsE)
