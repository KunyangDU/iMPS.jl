using TensorKit,JLD2,FiniteLattices

include("../../../src/iMPS.jl")
include("../model.jl")

Lx = 8
Ly = 1

Jx = -1
Jy = Jx
Jz = Jx

d = 2
D_MPS = 2^3

TruncErr = 1e-2
MaxIter = 20

Latt = load("examples/Heisenberg/data/$(Lx)x$(Ly)/Latt_$(Lx)x$(Ly).jld2")["Latt"]
lsE = load("examples/Heisenberg/data/$(Lx)x$(Ly)/lsE_D=$(D_MPS)_$(Lx)x$(Ly)_J=($((Jx,Jy,Jz))).jld2")["lsE"]
ψ = load("examples/Heisenberg/data/$(Lx)x$(Ly)/ψ_D=$(D_MPS)_$(Lx)x$(Ly)_J=($((Jx,Jy,Jz))).jld2")["ψ"]
H = load("examples/Heisenberg/data/$(Lx)x$(Ly)/H_$(Lx)x$(Ly)_J=($((Jx,Jy,Jz))).jld2")["H"]


σx = [0 1;1 0]
σy = [0 -1im;1im 0]
σz = [1 0;0 -1]
σs = [σx,σy, σz]

Sx = UnivMPO(size(Latt),σx)
Sy = UnivMPO(size(Latt),σy)
Sz = UnivMPO(size(Latt),σz)

Si = zeros(3,size(Latt))

for i in 1:size(Latt)
    for (ind,σi) in enumerate(σs)
        Szi = LocalMPO(size(Latt),σi,i)
        Si[ind,i] = QuantUniv(ψ,Szi)
    end
end

#QuantUniv(ψ,Sx),QuantUniv(ψ,Sy),QuantUniv(ψ,Sz)

Si

