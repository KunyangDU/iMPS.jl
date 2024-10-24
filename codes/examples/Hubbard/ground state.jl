include("../../src/iMPS.jl")
include("model.jl")


Lx = 8
Ly = 1
Latt = YCSqua(Lx,Ly)
t = 1
U = 8
d = 4
D_MPS = 2^5
lsμ = (U/4 - 4):0.5:(3*U/4 + 4)
@save "examples/Hubbard/data/lsμ_$(Lx)x$(Ly)_D_MPS=$(D_MPS)_t=$(t)_U=$(U).jld2" lsμ


LanczosLevel = 15
Nsweep = 3

for μ in 0.0
    @show μ
    H,D_MPO = compress(canonicalize(Hamiltonian(Latt;t=t,μ=μ,U=U)))
    ψ = RandMPS(Lx*Ly;d=d)
    ψ,lsE = sweepDMRG2(ψ,H,Nsweep,LanczosLevel,D_MPS)

    showQuantSweep(lsE;name="Eg sweep")
    @save "examples/Hubbard/data/ψ_$(Lx)x$(Ly)_D_MPS=$(D_MPS)_t=$(t)_μ=$(μ)_U=$(U).jld2" ψ
    @save "examples/Hubbard/data/lsE_$(Lx)x$(Ly)_D_MPS=$(D_MPS)_t=$(t)_μ=$(μ)_U=$(U).jld2" lsE

    ObsDict = let Obsf = ObserableForest()
        LocalSpace = Spin2Fermion
        for i in 1:size(Latt)
            addObs!(Obsf,LocalSpace.nup,i,"nup",nothing)
            addObs!(Obsf,LocalSpace.ndown,i,"ndown",nothing)
        end
    
        calObs(ψ, Obsf)
    end
    
    @save "examples/Hubbard/data/ObsDict_$(Lx)x$(Ly)_D_MPS=$(D_MPS)_t=$(t)_μ=$(μ)_U=$(U).jld2" ObsDict

end



