using TensorKit,JLD2,BenchmarkTools

include("../../../../src/iMPS.jl")
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
t = 2.0*abs(J)
Nt = 20
LanczosLevel = 20

lsψ,lst = sweepTDVP2(ψ,H,t,Nt,D_MPS,LanczosLevel;TruncErr = D_MPS)

@save "examples/MANUAL/Transverse Ising/data/tdvp/Impur_lsψ_D=$(D_MPS)_L=$(L)_J=$(J)_h=$(h).jld2" lsψ
@save "examples/MANUAL/Transverse Ising/data/tdvp/Impur_lst_D=$(D_MPS)_L=$(L)_J=$(J)_h=$(h).jld2" lst

#@benchmark sweepTDVP2(ψ,H,t,Nt,D_MPS;TruncErr = D_MPS)

