using TensorKit,JLD2,FiniteLattices

include("../../../src/iMPS.jl")
include("../model.jl")

lsLx = 4:2:20
Ly = 1

Jx = -1
Jy = Jx
Jz = Jx

d = 2
D_MPS = 2^5


for Lx in lsLx

    Latt = YCSqua(Lx,Ly)
    H = HamMPO(size(Latt);Jx=Jx,Jy=Jy,Jz=Jz,hz=0)
    ψ = RandMPS(size(Latt))

    LanczosLevel = 20
    Nsweep = 5

    ψ,lsE = sweepDMRG2(ψ,H,Nsweep,LanczosLevel,D_MPS)
    showQuantSweep(lsE;name = "Eg sweep")
    @save "examples/Heisenberg/data/$(Lx)x$(Ly)/Latt_$(Lx)x$(Ly).jld2" Latt
    @save "examples/Heisenberg/data/$(Lx)x$(Ly)/H_$(Lx)x$(Ly)_J=($((Jx,Jy,Jz))).jld2" H
    @save "examples/Heisenberg/data/$(Lx)x$(Ly)/ψ_D=$(D_MPS)_$(Lx)x$(Ly)_J=($((Jx,Jy,Jz))).jld2" ψ
    @save "examples/Heisenberg/data/$(Lx)x$(Ly)/lsE_D=$(D_MPS)_$(Lx)x$(Ly)_J=($((Jx,Jy,Jz))).jld2" lsE
end
# antiferro, ferro
# @show 0.4431*4*Jx*size(Latt)
