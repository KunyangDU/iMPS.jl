using TensorKit,JLD2,LinearAlgebra,FiniteLattices
include("../model.jl")
include("../../../src/iMPS.jl")

Lx = 7
Ly = 4
Latt = YCSqua(Lx,Ly)
@save "examples/Free Fermion/data/$(Lx)x$(Ly)/Latt_$(Lx)x$(Ly).jld2" Latt

d = 2

lsμ = -4.0:0.2:4.0
#lsμ = range(-2.0,2.0,Lx+1)
@save "examples/Free Fermion/data/$(Lx)x$(Ly)/lsμ_$(Lx)x$(Ly).jld2" lsμ

t = 1

D_MPS = 2^5
maxd = FindMaxDist(neighbor(Latt))
D_MPO = d*(2*maxd + 2)

LanczosLevel = D_MPO*d
Nsweep = 3

for μ in lsμ
    H = HamMPO(Latt;μ=μ)
    ψ = RandMPS(Lx*Ly)
    ψ,lsE1 = sweepDMRG2(ψ,H,Nsweep,LanczosLevel,D_MPS)
    ψ,lsE2 = sweepDMRG1(ψ,H,Nsweep,LanczosLevel,D_MPS)

    showQuantSweep(lsE2;name="Eg sweep2")
    showQuantSweep(lsE1;name="Eg sweep1")
    lsE = vcat(lsE1,lsE2)

    @save "examples/Free Fermion/data/$(Lx)x$(Ly)/ψ_D=$(D_MPS)_$(Lx)x$(Ly)_t=$(t)_μ=$(μ).jld2" ψ
    @save "examples/Free Fermion/data/$(Lx)x$(Ly)/lsE_D=$(D_MPS)_$(Lx)x$(Ly)_t=$(t)_μ=$(μ).jld2" lsE
end



