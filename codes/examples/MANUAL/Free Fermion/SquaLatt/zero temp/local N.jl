include("../../../../src/iMPS.jl")
include("../../model.jl")


Lx = 4
Ly = 4

t = 1

D_MPS = 2^6

Latt = load("examples/Free Fermion/data/Latt_$(Lx)x$(Ly).jld2")["Latt"]
lsμ = load("examples/Free Fermion/data/lsμ_$(Lx)x$(Ly).jld2")["lsμ"]

localN = LocalNMPO()
localNμ = zeros(Float64,length(lsμ),size(Latt))

for (i,μ) in enumerate(lsμ)
    ψ = load("examples/Free Fermion/data/ψ_D=$(D_MPS)_$(Lx)x$(Ly)_t=$(t)_μ=$(μ).jld2")["ψ"]
    localNμ[i,:] = Quant1(ψ,localN,D_MPS)
end
@save "examples/Free Fermion/data/localNμ_$(Lx)x$(Ly).jld2" localNμ


