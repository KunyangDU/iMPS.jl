using CairoMakie,JLD2,TensorKit,LaTeXStrings
include("../../src/MPSanalysis.jl")
include("../model.jl")

Lx = 6
Ly = 1

t = 1

D_MPS = 2^5

lsμ = load("../codes/examples/Free Fermion/data/$(Lx)x$(Ly)/lsμ_$(Lx)x$(Ly).jld2")["lsμ"]
Nμ = load("../codes/examples/Free Fermion/data/$(Lx)x$(Ly)/Nμ_D=$(D_MPS)_$(Lx)x$(Ly).jld2")["Nμ"]

nμ = Nμ / (Lx*Ly)
dμ = lsμ[2]-lsμ[1]


ind = 1:4:length(lsμ)
χ = diff(nμ[ind]) ./ diff(lsμ[ind])
centerμ = centralize(lsμ[ind])

theoμ = range(extrema(lsμ)...,100)
theonμ = @. (asin(theoμ/(2t))-asin(theoμ[1]/(2t))) / pi
theoχ = @. (1/pi) / sqrt((2*t)^2 - theoμ^2)

width,height = 0.9 .* (300,200)

fig = Figure()

axμ = Axis(fig[1,1],
xlabel = L"n",
ylabel = L"μ",
title = "$(Lx)x$(Ly) Squa Free Fermion: D=$(D_MPS)",
titlealign = :left,
width = width,height = height)
ylims!(axμ,1.1.*extrema(lsμ)...)
scatterlines!(axμ,nμ,lsμ)
lines!(axμ,theonμ,theoμ,color = :red,linewidth = 2.0)

axχ = Axis(fig[1,2],
xlabel = L"χ",
xticks = 0:0.1:1,
width = 0.2*width,height = height)
hideydecorations!(axχ,grid = false)

scatterlines!(axχ,χ,centerμ)
#lines!(axχ,theoχ,theoμ,color = :red)
ylims!(axχ,1.1.*extrema(lsμ)...)
xlims!(axχ,0,0.3)

resize_to_layout!(fig)

display(fig)

save("Free Fermion/figures/χ_D=$(D_MPS)_$(Lx)x$(Ly).pdf",fig)
