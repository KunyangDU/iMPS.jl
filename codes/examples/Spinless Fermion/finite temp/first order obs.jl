include("../../../src/iMPS.jl")
include("../model.jl")


Lx = 8
Ly = 1

t = 1
d = 2

D_MPO = 20

Latt = YCSqua(Lx,Ly)
#lsμ = load("examples/Free Fermion/data/finite temp/$(Lx)x$(Ly)/lsμ_$(Lx)x$(Ly)_D_MPO=$(D_MPO).jld2")["lsμ"]
for μ in [2.0]
lsβ = load("examples/Spinless Fermion/data/finite temp/lsβ_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_t=$(t)_μ=$(μ).jld2")["lsβ"]
lsρ = load("examples/Spinless Fermion/data/finite temp/lsρ_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_t=$(t)_μ=$(μ).jld2")["lsρ"]

H,D_MPOH = compress(canonicalize(Hamiltonian(Latt;t=t,μ=μ)))
N,D_MPON = compress(canonicalize(ParticleNumber(Latt)))

Fi = Vector{Float64}(undef,length(lsβ))
Ni = Vector{Float64}(undef,length(lsβ))
Ui = Vector{Float64}(undef,length(lsβ))
Si = Vector{Float64}(undef,length(lsβ))

lsβ1 = 2*lsβ
lsObsDict = load("examples/Spinless Fermion/data/lsObsDict_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_t=$(t)_μ=$(μ).jld2" )["lsObsDict"]

for (βi,β) in enumerate(lsβ)
    Z = Trace(lsρ[βi])
    Fi[βi] = - log(Z) / lsβ1[βi]
    Ni[βi] = Trace(lsρ[βi],N;Z=Z)
    Ui[βi] = Trace(lsρ[βi],H;Z=Z) + μ*Ni[βi]
    Si[βi] = lsβ1[βi]*(Ui[βi] - Fi[βi])
end

fi,ni,ui,si = map(x -> x / size(Latt),[Fi,Ni,Ui,Si])

@save "examples/Spinless Fermion/data/finite temp/lsβ1_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_t=$(t)_μ=$(μ).jld2" lsβ1
@save "examples/Spinless Fermion/data/finite temp/fi_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_t=$(t)_μ=$(μ).jld2" fi
@save "examples/Spinless Fermion/data/finite temp/ni_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_t=$(t)_μ=$(μ).jld2" ni
@save "examples/Spinless Fermion/data/finite temp/ui_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_t=$(t)_μ=$(μ).jld2" ui
@save "examples/Spinless Fermion/data/finite temp/si_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_t=$(t)_μ=$(μ).jld2" si

end



