using TensorKit,LinearAlgebra,JLD2

include("../../../src/iMPS.jl")
include("../model.jl")

L = 11

d = 2
D_MPS = 2^3

J = -1.0
h = -0.5
state = 6

H = HamMPO(L;J=J,h=h)
ψ = ImpurMPS(L,state)
#ψ = FerroMPS(L,"AFM")
t = 12.0*abs(J)
Nt = 60

lsψ,lst = sweepTDVP2(ψ,H,t,Nt,D_MPS;TruncErr = D_MPS)


@save "examples/Transverse Ising/data/tdvp/Impur_lsψ_D=$(D_MPS)_L=$(L)_J=$(J)_h=$(h).jld2" lsψ
@save "examples/Transverse Ising/data/tdvp/Impur_lst_D=$(D_MPS)_L=$(L)_J=$(J)_h=$(h).jld2" lst


