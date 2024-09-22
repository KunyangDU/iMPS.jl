include("../../../../../src/iMPS.jl")
include("../../model.jl")

Lx = 4
Ly = 4
Latt = YCSqua(Lx,Ly)
@save "examples/MANUAL/Free Fermion/data/$(Lx)x$(Ly)/Latt_$(Lx)x$(Ly).jld2" Latt

d = 2

lsμ = -2.0:0.2:2.0
@save "examples/MANUAL/Free Fermion/data/$(Lx)x$(Ly)/lsμ_$(Lx)x$(Ly).jld2" lsμ

t = 1

D_MPS = 2^4
maxd = FindMaxDist(neighbor(Latt))
D_MPO = d*(2*maxd + 2)

LanczosLevel = 20
Nsweep = 5
for μ in [1.0]
    @show μ
    H = HamMPO(Latt;μ=μ)
    ψ = RandMPS(Lx*Ly;d=d)
    ψ,lsE = sweepDMRG2(ψ,H,Nsweep,LanczosLevel,D_MPS)

    showQuantSweep(lsE;name="Eg sweep")
    
    @save "examples/MANUAL/Free Fermion/data/$(Lx)x$(Ly)/ψ_D=$(D_MPS)_$(Lx)x$(Ly)_t=$(t)_μ=$(μ).jld2" ψ
    @save "examples/MANUAL/Free Fermion/data/$(Lx)x$(Ly)/lsE_D=$(D_MPS)_$(Lx)x$(Ly)_t=$(t)_μ=$(μ).jld2" lsE
end
#= H = HamMPO(Latt;μ=μ)
ψ = RandMPS(Lx*Ly)
@benchmark sweepDMRG2(ψ,H,Nsweep,LanczosLevel,D_MPS) =#

