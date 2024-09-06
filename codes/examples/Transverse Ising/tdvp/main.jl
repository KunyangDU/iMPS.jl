using TensorKit,LinearAlgebra,JLD2

include("../../../src/iMPS.jl")
include("../model.jl")

L = 12

d = 2
D_MPS = 2^5

J = -1.0
h = 0.5
state = 6

H = IsingHam(L;J=J,h=h)
ψ = IsingLocalMPS(L,state,D_MPS)
t = 5.0*abs(J)
Nt = 100

lsψ,lst = sweepTDVP1(ψ,H,t,Nt,D_MPS)

#@save "trans Ising/data/tdvp/lsψ_D=$(D_MPS)_L=$(L)_J=$(J)_h=$(h).jld2" lsψ
#@save "trans Ising/data/tdvp/lst_D=$(D_MPS)_L=$(L)_J=$(J)_h=$(h).jld2" lst


