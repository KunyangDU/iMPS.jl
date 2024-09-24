include("../../../src/iMPS.jl")
include("model.jl")
include("../../MANUAL/Transverse Ising/model.jl")

Lx = 8
Ly = 1
J = 1
h = 0
hz = 0.0
Latt = YCSqua(Lx,Ly)

Nsweep = 3
D_MPS = 2^3

H,D_MPO = compress(canonicalize(Hamiltonian(Latt;J=J,h=h,hz=hz)))
LanczosLevel = 16

ψ = InitialRandΨ(Latt)

ψ,lsE = sweepDMRG2(ψ,H,Nsweep,LanczosLevel,D_MPS)

showQuantSweep(lsE)

ObsDict = let Obsf = ObserableForest()
    LocalSpace = Spin2
    for i in 1:size(Latt)
        addObs!(Obsf,LocalSpace.Sz,i,"Sz",nothing)
    end

    for pair in neighbor(Latt)
        addObs!(Obsf,LocalSpace.SzSz,pair,("Sz","Sz"),nothing)
    end

    for pair in neighbor(Latt;level = 2)
        addObs!(Obsf,LocalSpace.SzSz,pair,("Sz","Sz"),nothing;ObsName = "SzSzNN")
    end

    calObs(ψ, Obsf)
end

#= @save "examples/AUTO/Transverse Ising/data/ψ_$(Lx)x$(Ly)_D_MPS=$(D_MPS)_J=$(J).jld2" ψ
@save "examples/AUTO/Transverse Ising/data/lsE_$(Lx)x$(Ly)_D_MPS=$(D_MPS)_J=$(J).jld2" lsE =#

