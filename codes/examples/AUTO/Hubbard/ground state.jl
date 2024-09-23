include("../../../src/iMPS.jl")
include("model.jl")


Lx = 6
Ly = 1
Latt = YCSqua(Lx,Ly)
t = 1
U = 8
d = 4

lsμ = (U/4 - 4):0.5:(3*U/4 + 4)


D_MPS = 20
maxd = FindMaxDist(neighbor(Latt))
D_MPO = 4*maxd + 2

LanczosLevel = D_MPO*d
Nsweep = 3

for μ in lsμ
    @show μ
    H,D_MPO = compress(canonicalize(Hamiltonian(Latt;t=t,μ=μ,U=U)))
    ψ = RandMPS(Lx*Ly;d=d)
    ψ,lsE = sweepDMRG2(ψ,H,Nsweep,LanczosLevel,D_MPS)

    showQuantSweep(lsE;name="Eg sweep")

end

#= μ = 0.7
H = HubbardHam(Latt;t=t,μ=μ,U=U)

showBlockMPO(H[1]) =#
