include("../../../src/iMPS.jl")
include("model.jl")

Lx = 8
Ly = 1
Latt = YCSqua(Lx,Ly)
t = 1
μ = 0.0

lsμ = -2.0:0.2:2.0

Nsweep = 5
LanczosLevel = 15
D_MPS = 2^3

H,D_MPO = compress(canonicalize(Hamiltonian(Latt;t=t,μ=μ)))

ψ = InitialRandΨ(Latt)

ψ,lsE = sweepDMRG2(ψ,H,Nsweep,LanczosLevel,D_MPS)

showQuantSweep(lsE)

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

@save "examples/AUTO/Spinless Fermion/data/ψ_$(Lx)x$(Ly)_D_MPS=$(D_MPS)_t=$(t)_μ=$(μ).jld2" ψ
@save "examples/AUTO/Spinless Fermion/data/lsE_$(Lx)x$(Ly)_D_MPS=$(D_MPS)_t=$(t)_μ=$(μ).jld2" lsE
@save "examples/AUTO/Spinless Fermion/data/ObsDict_$(Lx)x$(Ly)_D_MPS=$(D_MPS)_t=$(t)_μ=$(μ).jld2" ObsDict
