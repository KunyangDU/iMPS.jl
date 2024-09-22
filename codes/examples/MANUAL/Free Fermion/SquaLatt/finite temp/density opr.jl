using BenchmarkTools
include("../../../../../src/iMPS.jl")
include("../../model.jl")


Lx = 12
Ly = 1

t = 1
lsμ = -2.0:0.4:2.0
d = 2

D_MPO = 32

LanczosLevel = 15

Latt = YCSqua(Lx,Ly)
#@save "examples/Free Fermion/data/finite temp/$(Lx)x$(Ly)/Latt_$(Lx)x$(Ly).jld2" Latt
#@save "examples/Free Fermion/data/finite temp/$(Lx)x$(Ly)/lsμ_$(Lx)x$(Ly)_D_MPO=$(D_MPO).jld2" lsμ

lsβ = vcat((3/2) .^ (-20:1:-1), 1:1/3:10)
β₀ = lsβ[1]
for μ in [1.0]
H = canonicalize(HamMPO(Latt;t=t,μ=μ,d=d))
ρ = SETTN(H,β₀,5,D_MPO)

lsρ = sweepTanTRG2(ρ,H,lsβ,D_MPO,LanczosLevel)

#@save "examples/Free Fermion/data/finite temp/$(Lx)x$(Ly)/lsβ_$(Lx)x$(Ly).jld2" lsβ
#@save "examples/Free Fermion/data/finite temp/$(Lx)x$(Ly)/lsρ_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_μ=$(μ).jld2" lsρ
end


