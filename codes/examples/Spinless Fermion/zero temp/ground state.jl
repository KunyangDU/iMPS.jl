include("../../../src/iMPS.jl")
include("../model.jl")

Lx = 16
Ly = 1
Latt = YCSqua(Lx,Ly)
t = 1
μ = 0.0
D_MPS = 2^5

lsμ = -3.0:0.2:3.0
@save "examples/Spinless Fermion/data/zero temp/lsμ_$(Lx)x$(Ly)_D_MPS=$(D_MPS)_t=$(t).jld2" lsμ

Nsweep = 3
LanczosLevel = 15

for μ in lsμ
H,D_MPO = compress(canonicalize(Hamiltonian(Latt;t=t,μ=μ)))

ψ = InitialRandΨ(Latt)

ψ,lsE = sweepDMRG2(ψ,H,Nsweep,LanczosLevel,D_MPS)

showQuantSweep(lsE)

@save "examples/Spinless Fermion/data/zero temp/ψ_$(Lx)x$(Ly)_D_MPS=$(D_MPS)_t=$(t)_μ=$(μ).jld2" ψ
@save "examples/Spinless Fermion/data/zero temp/lsE_$(Lx)x$(Ly)_D_MPS=$(D_MPS)_t=$(t)_μ=$(μ).jld2" lsE

ObsDict = let Obsf = ObserableForest()
    LocalSpace = SpinlessFermion
    for i in 1:size(Latt)
        addObs!(Obsf,LocalSpace.n,i,"n",nothing)
    end

    for pair in neighbor(Latt)
        addObs!(Obsf,LocalSpace.nn,pair,("n","n"),nothing)
    end

    calObs(ψ, Obsf)
end

@save "examples/Spinless Fermion/data/zero temp/ObsDict_$(Lx)x$(Ly)_D_MPS=$(D_MPS)_t=$(t)_μ=$(μ).jld2" ObsDict
end