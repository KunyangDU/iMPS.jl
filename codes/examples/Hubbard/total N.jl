

using TensorKit,JLD2,FiniteLattices
include("../../src/iMPS.jl")
include("../automodel.jl")

Lx = 6
Ly = 1

t = 1
U = 0
d = 4
μ = 0.0

D_MPS = 20

Latt = load("examples/Hubbard/data/$(Lx)x$(Ly)/Latt_$(Lx)x$(Ly).jld2")["Latt"]
lsμ = load("examples/Hubbard/data/$(Lx)x$(Ly)/lsμ_$(Lx)x$(Ly)_U=$(U).jld2")["lsμ"]
ψ = load("examples/Hubbard/data/$(Lx)x$(Ly)/ψ_D=$(D_MPS)_$(Lx)x$(Ly)_t=$(t)_U=$(U)_μ=$(μ).jld2")["ψ"]


#calObs(ψ,ParticleNumber(Latt))
ParticleNumber(Latt,1)
HubbardHam(Latt;t=1,U=8,returntree = true)


