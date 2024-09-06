using TensorKit,LinearAlgebra,JLD2

include("../../iMPS/iMPS.jl")
include("../model.jl")




L = 13

d = 2
D_MPS = 2^4

J = -1.0
h = 0.2

lsψ = load("trans Ising/data/tdvp/lsψ_D=$(D_MPS)_L=$(L)_J=$(J)_h=$(h).jld2")["lsψ"]
lst = load("trans Ising/data/tdvp/lst_D=$(D_MPS)_L=$(L)_J=$(J)_h=$(h).jld2")["lst"]

C = Vector{AbstractTensorMap}(undef,L)
DiffMi = zeros(Float64,length(lst),L)

MiMPO = TensorMap([1 0;0 -1],(ℂ^2)'→ (ℂ^2)')

for (it,t) in enumerate(lst)
    DiffMi[it,:] = Quant1(lsψ[it],MiMPO,D_MPS)
end

@save "trans Ising/data/tdvp/Mi_D=$(D_MPS)_L=$(L)_J=$(J)_h=$(h).jld2" DiffMi





