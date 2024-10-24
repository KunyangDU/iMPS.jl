include("../../../src/iMPS.jl")
include("../model.jl")


Lx = 8
Ly = 1

t = 1
lsμ = -2.0:0.4:2.0
d = 2

D_MPO = 20

LanczosLevel = 15

Latt = YCSqua(Lx,Ly)
@save "examples/Spinless Fermion/data/finite temp/lsμ_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_t=$(t).jld2" lsμ
L = size(Latt)
for μ in [2.0]

    lsβ = vcat((3/2) .^ (-20:1:-1), 1:1/3:10)
    β₀ = lsβ[1]

    H,D_MPO0 = compress(canonicalize(Hamiltonian(Latt;t=t,μ=μ)))
    ρ = SETTN(H,β₀,5,D_MPO)

    lsρ = sweepTanTRG2(ρ,H,lsβ,D_MPO,LanczosLevel)
    @save "examples/Spinless Fermion/data/finite temp/lsβ_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_t=$(t)_μ=$(μ).jld2" lsβ
    @save "examples/Spinless Fermion/data/finite temp/lsρ_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_t=$(t)_μ=$(μ).jld2" lsρ
    lsObsDict = Vector{Dict}(undef,length(lsβ))
    for (i,tempρ) in enumerate(lsρ)

        ObsDict = let Obsf = ObserableForest()
            LocalSpace = SpinlessFermion
            for i in 1:size(Latt)
                addObs!(Obsf,LocalSpace.n,i,"n",nothing)
            end

            calObs(tempρ, Obsf)
        end

        lsObsDict[i] = ObsDict

    end
    
    @save "examples/Spinless Fermion/data/lsObsDict_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_t=$(t)_μ=$(μ).jld2" lsObsDict


end


