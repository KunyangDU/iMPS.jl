include("../../../../src/iMPS.jl")
include("../../model.jl")


Lx = 6
Ly = 1

t = 1
d = 2

D_MPO = 32

Latt = load("examples/Free Fermion/data/finite temp/$(Lx)x$(Ly)/Latt_$(Lx)x$(Ly).jld2")["Latt"]
lsβ = load("examples/Free Fermion/data/finite temp/$(Lx)x$(Ly)/lsβ_$(Lx)x$(Ly).jld2")["lsβ"]
lsμ = load("examples/Free Fermion/data/finite temp/$(Lx)x$(Ly)/lsμ_$(Lx)x$(Ly)_D_MPO=$(D_MPO).jld2")["lsμ"]
for μ in [1.0]
lsρ = load("examples/Free Fermion/data/finite temp/$(Lx)x$(Ly)/lsρ_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_μ=$(μ).jld2")["lsρ"]

H = canonicalize(HamMPO(Latt;t=t,μ=μ,d=d))
N = canonicalize(NMPO(size(Latt)))

Fi = Vector{Float64}(undef,length(lsβ))
Ni = Vector{Float64}(undef,length(lsβ))
Ui = Vector{Float64}(undef,length(lsβ))
Si = Vector{Float64}(undef,length(lsβ))

lsβ1 = 2*lsβ

for (βi,β) in enumerate(lsβ)
    Z = Trace(lsρ[βi])
    Fi[βi] = - log(Z) / lsβ1[βi]
    Ni[βi] = Trace(lsρ[βi],N;Z=Z)
    Ui[βi] = Trace(lsρ[βi],H;Z=Z) + μ*Ni[βi]
    Si[βi] = lsβ1[βi]*(Ui[βi] - Fi[βi])
end

fi,ni,ui,si = map(x -> x / size(Latt),[Fi,Ni,Ui,Si])


@save "examples/Free Fermion/data/finite temp/$(Lx)x$(Ly)/lsβ1_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_μ=$(μ).jld2" lsβ1
@save "examples/Free Fermion/data/finite temp/$(Lx)x$(Ly)/fi_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_μ=$(μ).jld2" fi
@save "examples/Free Fermion/data/finite temp/$(Lx)x$(Ly)/ni_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_μ=$(μ).jld2" ni
@save "examples/Free Fermion/data/finite temp/$(Lx)x$(Ly)/ui_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_μ=$(μ).jld2" ui
@save "examples/Free Fermion/data/finite temp/$(Lx)x$(Ly)/si_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_μ=$(μ).jld2" si
end
