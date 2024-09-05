using TensorKit,LinearAlgebra,JLD2

include("../iMPS/iMPS.jl")
include("model.jl")

L = 12

J = -1.0
lsh = collect(-2.0:0.2:2.0)

D_MPS = 2^4


lsEg = Vector{Float64}(undef,length(lsh))
lsMz = Vector{Float64}(undef,length(lsh))

for (ih,h) in enumerate(lsh)
    ψ = load("trans Ising/data/ψ_D=$(D_MPS)_L=$(L)_J=$(J)_h=$(h).jld2")["ψ"]

    MM = IsingMagmom(L)
    reduM = reduHam(ψ,MM,1)

    H = IsingHam(L;J=J,h=h)
    reduH = reduHam(ψ,H,1)

    M = @tensor ψ[1][1,2]*reduM[1,2,3,4]*ψ[1]'[3,4]
    E = @tensor ψ[1][1,2]*reduH[1,2,3,4]*ψ[1]'[3,4]

    if abs(imag(M)) > 1e-5 || abs(imag(E)) > 1e-5
        @error "Quantity not real"
    else
        lsMz[ih] = real(M)
        lsEg[ih] = real(E)
    end
end

@save "trans Ising/data/lsh_D=$(D_MPS)_L=$(L)_J=$(J).jld2" lsh
@save "trans Ising/data/lsMz_D=$(D_MPS)_L=$(L)_J=$(J).jld2" lsMz
@save "trans Ising/data/lsEg_D=$(D_MPS)_L=$(L)_J=$(J).jld2" lsEg
