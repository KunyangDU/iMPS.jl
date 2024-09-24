using JLD2,TensorKit
include("../../../../src/iMPS.jl")
include("../model.jl")
include("../../../AUTO/Transverse Ising/model.jl")


Lx = 8
Ly = 1
Latt = YCSqua(Lx,Ly)

d = 2
D_MPS = 2^3

J = 1

h = 0.0
    
H,D_MPO = compress(canonicalize(HamMPO(size(Latt);J=J,h=h)))
ψ = RandMPS(size(Latt))

LanczosLevel = 16
Nsweep = 3

ψ,lsE = sweepDMRG2(ψ,H,Nsweep,LanczosLevel,D_MPS)

showQuantSweep(lsE;name="Eg sweep")



Obsf = ObserableForest()
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


#= @save "examples/MANUAL/Transverse Ising/data/dmrg/ψ_D=$(D_MPS)_L=$(L)_J=$(J)_h=$(h).jld2" ψ
@save "examples/MANUAL/Transverse Ising/data/dmrg/lsE_D=$(D_MPS)_L=$(L)_J=$(J)_h=$(h).jld2" lsE =#
