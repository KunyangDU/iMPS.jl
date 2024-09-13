using TensorKit,JLD2,LinearAlgebra,FiniteLattices
include("../model.jl")
include("../../../src/iMPS.jl")


Lx = 8
Ly = 4

t = 1

D_MPS = 2^3

Latt = load("examples/Free Fermion/data/$(Lx)x$(Ly)/Latt_$(Lx)x$(Ly).jld2")["Latt"]
lsμ = load("examples/Free Fermion/data/$(Lx)x$(Ly)/lsμ_$(Lx)x$(Ly).jld2")["lsμ"]

N = NMPO(size(Latt))
Nμ = Vector{Float64}(undef,length(lsμ))

for (i,μ) in enumerate(lsμ)
    ψ = load("examples/Free Fermion/data/$(Lx)x$(Ly)/ψ_D=$(D_MPS)_$(Lx)x$(Ly)_t=$(t)_μ=$(μ).jld2")["ψ"]
    Nμ[i] = QuantUniv(ψ,N)
end
@save "examples/Free Fermion/data/$(Lx)x$(Ly)/Nμ_D=$(D_MPS)_$(Lx)x$(Ly).jld2" Nμ


