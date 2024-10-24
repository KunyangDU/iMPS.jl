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
=======
D_MPO = 20
>>>>>>> 8ea8417fd317c4adb4f58a9cd6b4c299e7c2f40e

Latt = YCSqua(Lx,Ly)
lsβ2 = load("../codes/examples/Spinless Fermion/data/finite temp/lsβ2_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_t=$(t)_μ=$(μ).jld2")["lsβ2"]
ce = load("../codes/examples/Spinless Fermion/data/finite temp/ce_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_t=$(t)_μ=$(μ).jld2")["ce"]

width,height = 0.6 .* (600,300)

fig = Figure()

axce = Axis(fig[1,1],
ylabel = L"C_e/N",
xlabel = L"T\ /\ t",
xscale=log10,
title = "$(Lx)x$(Ly) Squa Free Fermion: D=$(D_MPO),μ=$(μ)",
titlealign = :left,
xminorticksvisible = true,
xminorgridvisible = true,
xminorticks = IntervalsBetween(5),
width = width,height = height,
)

scatterlines!(axce,1 ./ lsβ2,ce)

theoβ = 10 .^ (range(log10.(extrema(lsβ2))...,200))
theoce = HeatCapacity1(Latt,theoβ,μ)

lines!(axce,1 ./ theoβ,theoce,color = :red)

resize_to_layout!(fig)
display(fig)


save("Free Fermion/figures/second order_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_t=$(t)_μ=$(μ).pdf",fig)

