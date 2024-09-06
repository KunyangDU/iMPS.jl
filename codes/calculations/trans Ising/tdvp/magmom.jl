using TensorKit,LinearAlgebra,JLD2

include("../../iMPS/iMPS.jl")
include("../model.jl")




L = 12

d = 2
D_MPS = 2^4

J = -1.0
h = 0.0

lsψ = load("trans Ising/data/tdvp/lsψ_D=$(D_MPS)_L=$(L)_J=$(J)_h=$(h).jld2")["lsψ"]
lst = load("trans Ising/data/tdvp/lst_D=$(D_MPS)_L=$(L)_J=$(J)_h=$(h).jld2")["lst"]

lsC = Vector{AbstractTensorMap}(undef,L)
DiffMi = zeros(Float64,length(lst),L)


I = diagm(ones(2))
I0 = zeros(2,2)
σx = [0 1;1 0]
σz = [1 0;0 -1]
phys = ℂ^d


for (it,t) in enumerate(lst)
    ψ = lsψ[it]
    lsC[1] = ψ[1]
    for i in 1:L
        C = lsC[i]
        MiMPO = TensorMap(exp(-1im*h*σx*t)'*σz*exp(-1im*h*σx*t) .- σz,phys'→ phys')

        if i ==1
            Mi = @tensor C[1,2]*MiMPO[1,3]*C'[3,2]
        elseif i == L 
            Mi = @tensor C[1,2]*MiMPO[2,4]*C'[1,4]
        else
            Mi = @tensor C[1,2,3]*MiMPO[2,4]*C'[1,4,3]
        end

        DiffMi[it,i] = Mi
        if i < L
            lsC[i:i+1] = RightMove(ψ[i+1],lsC[i],D_MPS)
        end
        #DiffMiMPO = IsingLocalizedMagmom(L,i;t=t,h=h) .- IsingLocalizedMagmom(L,i;t=t,h=0)
    end
end


DiffMi





