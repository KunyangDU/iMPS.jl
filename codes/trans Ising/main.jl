using TensorKit,LinearAlgebra,JLD2

include("../iMPS/iMPS.jl")
include("model.jl")

L = 12

d = 2
D_MPO = 3
D_MPS = 128

lsParams = [(J = 0.0,h = -0.5),(J = 1.0,h = 0.)]

for params in lsParams
    
    H = tranIsingH(L,d,D_MPO;params...)
    ψ = initialMPS(L,d,D_MPS)
    p = ψ[:]

    LanczosLevel = 8
    Nsweep = 10

    ψ,lsE = sweep1(ψ,H,Nsweep,LanczosLevel,D_MPS)

    for i in eachindex(lsE)
        println(lsE[i])
    end

    @save "trans Ising/data/ψ_D=$(D_MPS)_L=$(L)_J=$(params[1])_h=$(params[2]).jld2" ψ
    @save "trans Ising/data/lsE_D=$(D_MPS)_L=$(L)_J=$(params[1])_h=$(params[2]).jld2" lsE
end


