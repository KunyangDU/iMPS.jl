using CairoMakie,JLD2,TensorKit,LaTeXStrings,FiniteLattices
include("../../src/MPSanalysis.jl")
include("../model.jl")

Lx = 8
Ly = 1

t = 1
μ = 2.0
d = 2

<<<<<<< HEAD
D_MPO = 32
Latt = load("../codes/examples/Free Fermion/data/finite temp/$(Lx)x$(Ly)/Latt_$(Lx)x$(Ly).jld2")["Latt"]
lsβ1 = load("../codes/examples/Free Fermion/data/finite temp/$(Lx)x$(Ly)/lsβ1_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_μ=$(μ).jld2")["lsβ1"]
lsβ2 = load("../codes/examples/Free Fermion/data/finite temp/$(Lx)x$(Ly)/lsβ2_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_μ=$(μ).jld2")["lsβ2"]
ni = load("../codes/examples/Free Fermion/data/finite temp/$(Lx)x$(Ly)/ni_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_μ=$(μ).jld2")["ni"]
fi = load("../codes/examples/Free Fermion/data/finite temp/$(Lx)x$(Ly)/fi_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_μ=$(μ).jld2")["fi"]
ui = load("../codes/examples/Free Fermion/data/finite temp/$(Lx)x$(Ly)/ui_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_μ=$(μ).jld2")["ui"]
si = load("../codes/examples/Free Fermion/data/finite temp/$(Lx)x$(Ly)/si_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_μ=$(μ).jld2")["si"]
=======
D_MPO = 20

Latt = YCSqua(Lx,Ly)

lsβ1 = load("../codes/examples/Spinless Fermion/data/finite temp/lsβ1_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_t=$(t)_μ=$(μ).jld2")["lsβ1"]
ni = load("../codes/examples/Spinless Fermion/data/finite temp/ni_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_t=$(t)_μ=$(μ).jld2")["ni"]
fi = load("../codes/examples/Spinless Fermion/data/finite temp/fi_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_t=$(t)_μ=$(μ).jld2")["fi"]
ui = load("../codes/examples/Spinless Fermion/data/finite temp/ui_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_t=$(t)_μ=$(μ).jld2")["ui"]
si = load("../codes/examples/Spinless Fermion/data/finite temp/si_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_t=$(t)_μ=$(μ).jld2")["si"]

lsβ2 = load("../codes/examples/Spinless Fermion/data/finite temp/lsβ2_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_t=$(t)_μ=$(μ).jld2")["lsβ2"]
>>>>>>> 8ea8417fd317c4adb4f58a9cd6b4c299e7c2f40e

width,height = 0.7 .* (600,150)

fig = Figure()

axN = Axis(fig[1,1],
ylabel = L"n",
xscale=log10,
#yticks = 0:0.25:1,
title = "$(Lx)x$(Ly) Squa Free Fermion: D=$(D_MPO),μ=$(μ)",
titlealign = :left,
xminorticksvisible = true,
xminorgridvisible = true,
xminorticks = IntervalsBetween(5),
width = width,height = height,
)

scatterlines!(axN,1 ./ lsβ1,ni)
#ylims!(axN,0.4,0.6)

axU = Axis(fig[2,1],
ylabel = L"u",
xscale=log10,
xminorticksvisible = true,
xminorgridvisible = true,
xminorticks = IntervalsBetween(5),
width = width,height = height,
)

scatterlines!(axU,1 ./ lsβ1,ui)

axF = Axis(fig[3,1],
ylabel = L"f",
xlabel = L"T\ /\ t",
xscale=log10,
xminorticksvisible = true,
xminorgridvisible = true,
xminorticks = IntervalsBetween(5),
width = width,height = height,
)
scatterlines!(axF,1 ./ lsβ1,fi)


theoβ = 10 .^ (range(log10.(extrema(lsβ1))...,300))
#theoβ = 10 .^ (-8.0:1:-4.0)
theof = FreeEnergy1(Latt,theoβ,μ)
theou = InternalEnergy1(Latt,theoβ,μ)
theon = ParticleNumber1(Latt,theoβ,μ)
lines!(axF,1 ./ theoβ,theof,color = :red)
lines!(axU,1 ./ theoβ,theou,color = :red)
lines!(axN,1 ./ theoβ,theon,color = :red)
#ylims!(axF,-1.5,-0.4)


hidexdecorations!(axN,grid = false)
hidexdecorations!(axU,grid = false)

resize_to_layout!(fig)
display(fig)

save("Free Fermion/figures/fist order_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_t=$(t)_μ=$(μ).pdf",fig)

