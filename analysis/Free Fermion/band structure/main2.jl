using CairoMakie,JLD2,TensorKit,LaTeXStrings
include("../../src/MPSanalysis.jl")
include("../model.jl")

Lx = 8
Ly = 4

t = 1

D_MPS = 2^3

lsμ = load("../codes/examples/Free Fermion/data/$(Lx)x$(Ly)/lsμ_$(Lx)x$(Ly).jld2")["lsμ"]
χk = load("../codes/examples/Free Fermion/data/$(Lx)x$(Ly)/χk_D=$(D_MPS)_$(Lx)x$(Ly).jld2")["χk"]
kr = load("../codes/examples/Free Fermion/data/$(Lx)x$(Ly)/χkr_D=$(D_MPS)_$(Lx)x$(Ly).jld2")["kr"]

centerμ = centralize(lsμ)

width,height = 0.9 .* (300,200)


xtickvalues = [0,pi,2*pi,(2+sqrt(2))*pi]
xticklabels = [L"(0,0)",L"(\pi,0)",L"(\pi,\pi)",L"(0,0)"]

fig = Figure()
ax = Axis(fig[1,1],
ylabel = L"\mu /t",
xlabel = L"k_r",
title = "$(Lx)x$(Ly) Squa Free Fermion: D=$(D_MPS)",
xticks = (xtickvalues,xticklabels),
width = width,height = height)
hm = heatmap!(ax,kr,centerμ,χk,colormap = :hot)
Colorbar(fig[1,2],hm,
label = L"\chi_k = \partial n_k / \partial \mu")

ipath = [0 pi pi 0;0 0 pi 0]
kvecpath = vrange(ipath;eachstep = 20)
theokr = pathlength(kvecpath)
band = [SquaBand(kv) for kv in collect.(eachcol(kvecpath))]
lines!(ax,theokr,band,color = :black)

#= totalNk = load("../codes/examples/Free Fermion/data/$(Lx)x$(Ly)/totalNk_D=$(D_MPS)_$(Lx)x$(Ly).jld2")["totalNk"]
axmu = Axis(fig[1,3])
scatterlines!(axmu,lsμ,totalNk)
scatterlines!(axmu,lsμ,Nμ)
 =#

resize_to_layout!(fig)
display(fig)

save("Free Fermion/figures/bandstructure_D=$(D_MPS)_$(Lx)x$(Ly).pdf",fig)

totalNk[end]




