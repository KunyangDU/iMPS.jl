include("../../../../src/iMPS.jl")
include("../../model.jl")


Lx = 8
Ly = 1

t = 1
d = 2

D_MPO = 2^6

Latt = load("examples/Free Fermion/data/finite temp/$(Lx)x$(Ly)/Latt_$(Lx)x$(Ly).jld2")["Latt"]
lsβ = load("examples/Free Fermion/data/finite temp/$(Lx)x$(Ly)/lsβ_$(Lx)x$(Ly).jld2")["lsβ"]
lsμ = load("examples/Free Fermion/data/finite temp/$(Lx)x$(Ly)/lsμ_$(Lx)x$(Ly)_D_MPO=$(D_MPO).jld2")["lsμ"]
for μ in [0.0]
lsρ = load("examples/Free Fermion/data/finite temp/$(Lx)x$(Ly)/lsρ_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_μ=$(μ).jld2")["lsρ"]
H = canonicalize(HamMPO(Latt;t=t,μ=μ,d=d))
N = canonicalize(NMPO(size(Latt)))
Ni = Vector{Float64}(undef,length(lsβ))
Fi = Vector{Float64}(undef,length(lsβ))
Ui = Vector{Float64}(undef,length(lsβ))
for (βi,β) in enumerate(lsβ)
    @assert Trace(lsρ[βi]) ≈ 1
    Fi[βi] = Trace(lsρ[βi],H)
    Ni[βi] = Trace(lsρ[βi],N)
end
Ui = Fi + μ*Ni

ni,fi,ui = map(x -> x / size(Latt),[Ni,Fi,Ui])

lsβ1 = 2*lsβ

@save "examples/Free Fermion/data/finite temp/$(Lx)x$(Ly)/lsβ1_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_μ=$(μ).jld2" lsβ1
@save "examples/Free Fermion/data/finite temp/$(Lx)x$(Ly)/ni_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_μ=$(μ).jld2" ni
@save "examples/Free Fermion/data/finite temp/$(Lx)x$(Ly)/fi_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_μ=$(μ).jld2" fi
@save "examples/Free Fermion/data/finite temp/$(Lx)x$(Ly)/ui_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_μ=$(μ).jld2" ui
end
