using TensorKit,JLD2,LinearAlgebra,FiniteLattices
include("../model.jl")
include("../../../../src/iMPS.jl")

# 赋值问题，统一为[A,B],err形式
# truncerr问题，为什么TDVP后面不再涨？

L = 12

d = 2
D_MPO = 3
D_MPS = 2^5

J = 1.0

H = HamMPO(L;J=J)

ψ = RandMPS(L)
LanczosLevel = 16
Nsweep = 3

ψ,lsE1 = sweepDMRG2(ψ,H,Nsweep,LanczosLevel,D_MPS)
k = [1,0]
CK = 



Hψ = VariContract(H,ψ,D_MPS;MaxIter = 3)
InnerProd(ψ,Hψ)


