include("../../src/iMPS.jl")
include("model.jl")
Lx = 21
Ly = 1
J = 1
h = 0.3
hz = 0.0
Latt = YCSqua(Lx,Ly)

Nsweep = 3
D_MPS = 2^5

H,D_MPO = compress(canonicalize(Hamiltonian(Latt;J=J,h=h,hz=hz)))
LanczosLevel = 16
name = "Center"
ψ = ManualMPS(Latt,name)

t = 40*J
Nt = 200

lsψ,lst = TDVP2!(ψ,H,t,Nt,D_MPS,LanczosLevel;TruncErr = D_MPS)

lsObsDict = Vector{Dict}(undef,length(lst))

for i in eachindex(lst)
    ObsDict = let Obsf = ObserableForest()
        LocalSpace = Spin2
        for i in 1:size(Latt)
            addObs!(Obsf,LocalSpace.Sz,i,"Sz",nothing)
        end
    
        calObs(lsψ[i], Obsf)
    end
    
    lsObsDict[i] = ObsDict
end

@save "examples/Transverse Ising/data/lsψ_$(Lx)x$(Ly)_D_MPS=$(D_MPS)_J=$(J)_h=$(h)_time=$(t)_name=$(name).jld2" lsψ
@save "examples/Transverse Ising/data/lst_$(Lx)x$(Ly)_D_MPS=$(D_MPS)_J=$(J)_h=$(h)_time=$(t)_name=$(name).jld2" lst
@save "examples/Transverse Ising/data/lsObsDict_$(Lx)x$(Ly)_D_MPS=$(D_MPS)_J=$(J)_h=$(h)_time=$(t)_name=$(name).jld2" lsObsDict
