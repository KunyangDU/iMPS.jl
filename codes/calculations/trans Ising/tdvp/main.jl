using TensorKit,LinearAlgebra,JLD2

include("../../iMPS/iMPS.jl")
include("../model.jl")

L = 12

d = 2
D_MPS = 2^4

J = -1.0
h = 0.0

H = IsingHam(L;J=J,h=h)
ψ = RandMPS(L,d,D_MPS)

t = 5.0
Nt = 50

lsψ,lst = sweepTDVP1(ψ,H,t,Nt,D_MPS)

@save "trans Ising/data/tdvp/lsψ_D=$(D_MPS)_L=$(L)_J=$(J)_h=$(h).jld2" lsψ
@save "trans Ising/data/tdvp/lst_D=$(D_MPS)_L=$(L)_J=$(J)_h=$(h).jld2" lst


