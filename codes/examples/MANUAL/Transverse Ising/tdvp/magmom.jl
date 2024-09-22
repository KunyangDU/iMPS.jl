using TensorKit,JLD2

include("../../../src/iMPS.jl")
include("../model.jl")


L = 11

d = 2
D_MPS = 2^3

J = -1.0
h = -0.5

lsψ = load("examples/Transverse Ising/data/tdvp/Impur_lsψ_D=$(D_MPS)_L=$(L)_J=$(J)_h=$(h).jld2")["lsψ"]
lst = load("examples/Transverse Ising/data/tdvp/Impur_lst_D=$(D_MPS)_L=$(L)_J=$(J)_h=$(h).jld2")["lst"]

Mi = zeros(Float64,length(lst),L)

MiMPO = LocalMagmomMPO()
for (it,t) in enumerate(lst)
    Mi[it,:] = Quant1(lsψ[it],MiMPO,D_MPS)
end

@save "examples/Transverse Ising/data/tdvp/Impur_Mi_D=$(D_MPS)_L=$(L)_J=$(J)_h=$(h).jld2" Mi





