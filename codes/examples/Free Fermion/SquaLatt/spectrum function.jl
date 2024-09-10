using TensorKit,JLD2,LinearAlgebra,FiniteLattices
include("../model.jl")
include("../../../src/iMPS.jl")

# 赋值问题，统一为[A,B],err形式
# truncerr问题

Lx = 6
Ly = 1
Latt = YCSqua(Lx,Ly)

d = 2

t = 1

D_MPS = 2^6
maxd = FindMaxDist(neighbor(Latt))
D_MPO = d*(2*maxd + 2)

LanczosLevel = D_MPO*d
Nsweep = 3

H = HamMPO(Latt)
ψ = RandMPS(Lx*Ly)
ψ,lsE = sweepDMRG2(ψ,H,Nsweep,LanczosLevel,D_MPS)

Gt,lst = GreenFuncTDVP2(ψ,H,0.1,1e-3,50,D_MPS)


