include("../../src/iMPS.jl")
include("model.jl")

Lx = 8
Ly = 1
J = -1
Latt = YCSqua(Lx,Ly)

Nsweep = 5
LanczosLevel = 15
D_MPS = 2^4

H,D_MPO = compress(canonicalize(Hamiltonian(Latt;Jx=J,h=0)))
ψ = InitialRandΨ(Latt)

ψ,lsE = sweepDMRG2(ψ,H,Nsweep,LanczosLevel,D_MPS)

showQuantSweep(lsE)
ObsDict = let Obsf = ObserableForest()
    LocalSpace = Spin2
    for i in 1:size(Latt)
        addObs!(Obsf,LocalSpace.Sx,i,"Sx",nothing)
        addObs!(Obsf,LocalSpace.Sy,i,"Sy",nothing)
        addObs!(Obsf,LocalSpace.Sz,i,"Sz",nothing)
    end

    for pair in neighbor(Latt)
        addObs!(Obsf,LocalSpace.SxSx,pair,("Sx","Sx"),nothing)
        addObs!(Obsf,LocalSpace.SySy,pair,("Sy","Sy"),nothing)
        addObs!(Obsf,LocalSpace.SzSz,pair,("Sz","Sz"),nothing)
    end

    calObs(ψ, Obsf)
end

@save "examples/Heisenberg/data/ψ_$(Lx)x$(Ly)_D_MPS=$(D_MPS)_J=$(J).jld2" ψ
@save "examples/Heisenberg/data/lsE_$(Lx)x$(Ly)_D_MPS=$(D_MPS)_J=$(J).jld2" lsE
@save "examples/Heisenberg/data/ObsDict_$(Lx)x$(Ly)_D_MPS=$(D_MPS)_J=$(J).jld2" ObsDict

