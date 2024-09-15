using TensorKit,JLD2,FiniteLattices
include("../model.jl")
include("../../../src/iMPS.jl")

Lx = 6
Ly = 1
Latt = YCSqua(Lx,Ly)
@save "examples/Hubbard/data/$(Lx)x$(Ly)/Latt_$(Lx)x$(Ly).jld2" Latt
t = 1
U = 8
d = 4

lsμ = -4.0:0.2:4.0
@save "examples/Hubbard/data/$(Lx)x$(Ly)/lsμ_$(Lx)x$(Ly).jld2" lsμ


D_MPS = 20
maxd = FindMaxDist(neighbor(Latt))
D_MPO = 4*maxd + 2

LanczosLevel = D_MPO*d
Nsweep = 5

for μ in [0.0]
    @show μ
    H = HamMPO(Latt;t=t,μ=μ,U=U,d=d)
    ψ = RandMPS(Lx*Ly;d=d)
    ψ,lsE = sweepDMRG2(ψ,H,Nsweep,LanczosLevel,D_MPS)

    showQuantSweep(lsE;name="Eg sweep")

    @save "examples/Hubbard/data/$(Lx)x$(Ly)/ψ_D=$(D_MPS)_$(Lx)x$(Ly)_t=$(t)_U=$(U)_μ=$(μ).jld2" ψ
    @save "examples/Hubbard/data/$(Lx)x$(Ly)/lsE_D=$(D_MPS)_$(Lx)x$(Ly)_t=$(t)_U=$(U)_μ=$(μ).jld2" lsE
end



