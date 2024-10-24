using FiniteLattices
include("../../../src/iMPS.jl")
include("../model.jl")


Lx = 8
Ly = 1
Latt = YCSqua(Lx,Ly)
@save "examples/Hubbard/data/$(Lx)x$(Ly)/Latt_$(Lx)x$(Ly).jld2" Latt
t = 1
U = 0
d = 4

lsμ = (U/4 - 4):0.5:(3*U/4 + 4)
@save "examples/Hubbard/data/$(Lx)x$(Ly)/lsμ_$(Lx)x$(Ly)_U=$(U).jld2" lsμ


D_MPS = 20
maxd = FindMaxDist(neighbor(Latt))
D_MPO = 4*maxd + 2

LanczosLevel = D_MPO*d
Nsweep = 3

for μ in lsμ
    @show μ
    H = HamMPO(Latt;t=t,μ=μ,U=U,d=d)
    ψ = RandMPS(Lx*Ly;d=d)
    ψ,lsE = sweepDMRG2(ψ,H,Nsweep,LanczosLevel,D_MPS)

    showQuantSweep(lsE;name="Eg sweep")

    @save "examples/Hubbard/data/$(Lx)x$(Ly)/ψ_D=$(D_MPS)_$(Lx)x$(Ly)_t=$(t)_U=$(U)_μ=$(μ).jld2" ψ
    @save "examples/Hubbard/data/$(Lx)x$(Ly)/lsE_D=$(D_MPS)_$(Lx)x$(Ly)_t=$(t)_U=$(U)_μ=$(μ).jld2" lsE
end



