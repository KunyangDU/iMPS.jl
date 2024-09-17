using TensorKit,JLD2,FiniteLattices,CairoMakie

include("model.jl")
include("../../src/iMPS.jl")

Lx = 4
Ly = 1
Latt = YCSqua(Lx,Ly)
d = 4
μ= 2.0
t=1
U = 8

D_MPS = 20

ψ = load("examples/Hubbard/data/$(Lx)x$(Ly)/ψ_D=$(D_MPS)_$(Lx)x$(Ly)_t=$(t)_U=$(U)_μ=$(μ).jld2")["ψ"]
Latt = load("examples/Hubbard/data/$(Lx)x$(Ly)/Latt_$(Lx)x$(Ly).jld2")["Latt"]

kv = [pi,0]
oprs = [aup,F*adown]
opr = 
ckup = KOprMPO(Latt,oprs[1],kv,-1;d=d,string=F)
ckupdagg = KOprMPO(Latt,collect(oprs[1]'),kv,1;d=d,string=F)
ckdown = KOprMPO(Latt,oprs[2],kv,1;d=d,string=F)
ckdowndagg = KOprMPO(Latt,collect(oprs[2]'),kv,1;d=d,string=F)

Nup = UnivMPO(size(Latt),oprs[1]'*oprs[1])
Ndown = UnivMPO(size(Latt),oprs[2]'*oprs[2])

ckupψ = VariContract(ckup,ψ,D_MPS;d=d)
ckdownψ = VariContract(ckdown,ψ,D_MPS;d=d)

InnerProd(ckupψ,ckupψ),InnerProd(ckdownψ,ckdownψ),QuantUniv(ψ,Nup),QuantUniv(ψ,Ndown),QuantUniv(ψ,Nup)+QuantUniv(ψ,Ndown)


