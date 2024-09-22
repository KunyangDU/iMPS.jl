using .iMPS
using FiniteLattices
include("model.jl")

Lx = 6
Ly = 1
Latt = YCSqua(Lx,Ly)
t = 1
μ = 0.0

lsμ = -2.0:0.2:2.0

Nsweep = 5
LanczosLevel = 15
D_MPS = 2^3

H = Hamiltonian(Latt;t=t,μ=μ)

ψ = SpinlessFermionRandMPS(Latt)

ψ,lsE = sweepDMRG2(ψ,H,Nsweep,LanczosLevel,D_MPS)

showQuantSweep(lsE)

calObs(ψ,ParticleNumber(Latt))

