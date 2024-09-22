using TensorKit,JLD2,FiniteLattices,CairoMakie

include("model.jl")
include("../../src/iMPS.jl")

Lx = 6
Ly = 1
Latt = YCSqua(Lx,Ly)
d = 2
μ= 0.0
t=1

D_MPS = 2^3

ψ = load("examples/Free Fermion/data/$(Lx)x$(Ly)/ψ_D=$(D_MPS)_$(Lx)x$(Ly)_t=$(t)_μ=$(μ).jld2")["ψ"]
Latt = load("examples/Free Fermion/data/$(Lx)x$(Ly)/Latt_$(Lx)x$(Ly).jld2")["Latt"]

kv = [0,0]
ck = KOprMPO(Latt,a,kv,-1;d=d,string = F)
ckdagg = KOprMPO(Latt,a⁺,kv,1;d=d,string = F)

N = NMPO(size(Latt))

ckψ = VariContract(ck,ψ,D_MPS)
ckdaggψ = VariContract(ckdagg,ψ,D_MPS)

InnerProd(ckψ,ckψ),InnerProd(ckdaggψ,ckdaggψ),QuantUniv(ckψ,N),QuantUniv(ckdaggψ,N),QuantUniv(ψ,N)


