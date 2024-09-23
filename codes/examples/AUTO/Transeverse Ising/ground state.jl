include("../../../src/iMPS.jl")
include("model.jl")

Lx = 4
Ly = 1
J = 1
Latt = YCSqua(Lx,Ly)

Nsweep = 5
LanczosLevel = 15
D_MPS = 2^3

H,D_MPO = compress(canonicalize(Hamiltonian(Latt;J=J,h=0)))

ψ = InitialRandΨ(Latt)

ψ,lsE = sweepDMRG2(ψ,H,Nsweep,LanczosLevel,D_MPS)

showQuantSweep(lsE)
@save "examples/AUTO/Transeverse Ising/data/ψ_$(Lx)x$(Ly)_D_MPS=$(D_MPS)_J=$(J).jld2" ψ
@save "examples/AUTO/Transeverse Ising/data/lsE_$(Lx)x$(Ly)_D_MPS=$(D_MPS)_J=$(J).jld2" lsE

