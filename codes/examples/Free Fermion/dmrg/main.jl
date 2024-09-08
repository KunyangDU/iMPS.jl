using TensorKit,JLD2,LinearAlgebra,FiniteLattices
include("../model.jl")
include("../../../src/iMPS.jl")

Lx = 6
Ly = 2
Latt = YCSqua(Lx,Ly)
@save "examples/Free Fermion/data/Latt_$(Lx)x$(Ly).jld2" Latt

d = 2

lsμ = -2.0:0.2:2.0
t = 1

Nsweep = 3
D_MPS = 2^6
maxd = FindMaxDist(neighbor(Latt))
D_MPO = d*(2*maxd + 2)

LanczosLevel = D_MPO*d
Nsweep = 5


H = HamMPO(Latt;μ=μ)


for μ in lsμ
    ψ = RandMPS(Lx*Ly)
    ψ,lsE1 = sweepDMRG2(ψ,H,Nsweep,LanczosLevel,D_MPS)
    ψ,lsE2 = sweepDMRG1(ψ,H,Nsweep,LanczosLevel,D_MPS)

    showQuantSweep(lsE2;name="Eg sweep2")
    showQuantSweep(lsE1;name="Eg sweep1")
    lsE = vcat(lsE1,lsE2)

    @save "examples/Free Fermion/data/ψ_D=$(D_MPS)_$(Lx)x$(Ly)_t=$(t)_μ=$(μ).jld2" ψ
    @save "examples/Free Fermion/data/lsE_D=$(D_MPS)_$(Lx)x$(Ly)_t=$(t)_μ=$(μ).jld2" lsE
end
