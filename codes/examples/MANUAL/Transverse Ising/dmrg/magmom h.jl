using TensorKit,LinearAlgebra,JLD2

include("../../../src/iMPS.jl")
include("../model.jl")

L = 12

d = 2
D_MPO = 3
D_MPS = 2^5

J = 1.0
lsh = collect(-2.0:0.1:2.0)

lsEg = Vector{Float64}(undef,length(lsh))
lsMz = Vector{Float64}(undef,length(lsh))

MzMPO = MagmomMPO(L)


for (ih,h) in enumerate(lsh)
    ψ = load("examples/Transverse Ising/data/dmrg/ψ_D=$(D_MPS)_L=$(L)_J=$(J)_h=$(h).jld2")["ψ"]

    H = HamMPO(L;J=J,h=h)

    lsMz[ih] = QuantUniv(ψ,MzMPO)
    lsEg[ih] = QuantUniv(ψ,H)

end

@save "examples/Transverse Ising/data/dmrg/lsh_D=$(D_MPS)_L=$(L)_J=$(J).jld2" lsh
@save "examples/Transverse Ising/data/dmrg/lsMz_D=$(D_MPS)_L=$(L)_J=$(J).jld2" lsMz
@save "examples/Transverse Ising/data/dmrg/lsEg_D=$(D_MPS)_L=$(L)_J=$(J).jld2" lsEg
