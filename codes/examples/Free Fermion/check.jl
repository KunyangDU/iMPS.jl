using TensorKit,JLD2,LinearAlgebra,FiniteLattices,CairoMakie
include("model.jl")
include("../Transverse Ising/model.jl")

include("../../src/iMPS.jl")

# 赋值问题，统一为[A,B],err形式
# truncerr问题，为什么TDVP后面不再涨？

Lx = 4
Ly = 4
d = 2
μ = 3.
t=1

D_MPS = 2^3

Latt = load("examples/Free Fermion/data/$(Lx)x$(Ly)/Latt_$(Lx)x$(Ly).jld2")["Latt"]

ψ = load("examples/Free Fermion/data/$(Lx)x$(Ly)/ψ_D=$(D_MPS)_$(Lx)x$(Ly)_t=$(t)_μ=$(μ).jld2")["ψ"]
kv = [0,0]
ck = CKMPO(Latt,kv)
ckdagg = CKdaggMPO(Latt,kv)
ckψ = VariContract(ck,ψ,D_MPS)
ckdaggψ = VariContract(ckdagg,ψ,D_MPS)

nkψ = VariContract(ckdagg,ckψ,D_MPS)
N = NMPO(size(Latt))

QuantUniv(ψ,N),InnerProd(ckψ,ckψ),InnerProd(ψ,nkψ)



