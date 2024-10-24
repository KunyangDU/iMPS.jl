using CairoMakie,JLD2,TensorKit,LaTeXStrings,FiniteLattices
include("../../src/MPSanalysis.jl")
include("../model.jl")

Lx = 6
Ly = 1

t = 1
μ = 1.0
d = 2

D_MPO = 32

Latt = load("../codes/examples/Free Fermion/data/finite temp/$(Lx)x$(Ly)/Latt_$(Lx)x$(Ly).jld2")["Latt"]
lsβ2 = load("../codes/examples/Free Fermion/data/finite temp/$(Lx)x$(Ly)/lsβ2_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_μ=$(μ).jld2")["lsβ2"]
ce = load("../codes/examples/Free Fermion/data/finite temp/$(Lx)x$(Ly)/ce_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_μ=$(μ).jld2")["ce"]

width,height = 0.7 .* (600,450)

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


save("Free Fermion/figures/second order_$(Lx)x$(Ly)_D=$(D_MPO)_μ=$(μ).pdf",fig)

