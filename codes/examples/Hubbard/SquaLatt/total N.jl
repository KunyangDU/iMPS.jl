using TensorKit,JLD2,FiniteLattices
include("../../../src/iMPS.jl")
include("../model.jl")

Lx = 8
Ly = 1

t = 1
U = 8
d = 4

D_MPS = 20

Latt = load("examples/Hubbard/data/$(Lx)x$(Ly)/Latt_$(Lx)x$(Ly).jld2")["Latt"]
lsμ = load("examples/Hubbard/data/$(Lx)x$(Ly)/lsμ_$(Lx)x$(Ly)_U=$(U).jld2")["lsμ"]


oprs = [aup,F*adown]
Nup = UnivMPO(size(Latt),oprs[1]'*oprs[1])
Ndown = UnivMPO(size(Latt),oprs[2]'*oprs[2])

Nμ = zeros(length(lsμ),2)

for (i,μ) in enumerate(lsμ)
    for (Ni,No) in enumerate([Nup,Ndown])
        ψ = load("examples/Hubbard/data/$(Lx)x$(Ly)/ψ_D=$(D_MPS)_$(Lx)x$(Ly)_t=$(t)_U=$(U)_μ=$(μ).jld2")["ψ"]
        Nμ[i,Ni] = QuantUniv(ψ,No)
    end
end
@save "examples/Hubbard/data/$(Lx)x$(Ly)/Nμ_D=$(D_MPS)_$(Lx)x$(Ly)_U=$(U).jld2" Nμ


