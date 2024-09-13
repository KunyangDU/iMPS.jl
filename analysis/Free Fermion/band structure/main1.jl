using CairoMakie,JLD2,TensorKit,LaTeXStrings
include("../../src/MPSanalysis.jl")
include("../model.jl")

Lx = 6
Ly = 1

t = 1

D_MPS = 2^5

lsμ = load("../codes/examples/Free Fermion/data/$(Lx)x$(Ly)/lsμ_$(Lx)x$(Ly).jld2")["lsμ"]
χk = load("../codes/examples/Free Fermion/data/$(Lx)x$(Ly)/χk_D=$(D_MPS)_$(Lx)x$(Ly).jld2")["χk"]
kr = load("../codes/examples/Free Fermion/data/$(Lx)x$(Ly)/χkr_D=$(D_MPS)_$(Lx)x$(Ly).jld2")["kr"]
totalNk = load("../codes/examples/Free Fermion/data/$(Lx)x$(Ly)/totalNk_D=$(D_MPS)_$(Lx)x$(Ly).jld2")["totalNk"]
Nμ = load("../codes/examples/Free Fermion/data/$(Lx)x$(Ly)/Nμ_D=$(D_MPS)_$(Lx)x$(Ly).jld2")["Nμ"]


centerμ = centralize(lsμ)

xtickvalues = [0,pi/2,pi,3*pi/2,2*pi]
xticklabels = [L"-\pi",L"-\pi/2",L"0",L"\pi/2",L"\pi"]

fig = Figure()
ax = Axis(fig[1,1],
ylabel = L"\mu /t",
xlabel = L"k_r",
xticks = (xtickvalues,xticklabels),
title = "$(Lx)x$(Ly) Squa Free Fermion: D=$(D_MPS)",
width = width,height = height)
hm = heatmap!(ax,kr,centerμ,χk)
Colorbar(fig[1,2],hm,
label = L"\chi_k = \partial n_k / \partial \mu")

ipath = [-pi pi;0 0]
kvecpath = vrange(ipath;eachstep = 100) 
theokr = pathlength(kvecpath)
band = [SquaBand(kv)-2 for kv in collect.(eachcol(kvecpath))]
lines!(ax,theokr,band,color = :red,linestyle = :dash)
for xv in xtickvalues
    lines!(ax,[xv,xv],collect(extrema(lsμ)),color = :black)
end

resize_to_layout!(fig)
display(fig)

save("Free Fermion/figures/bandstructure_D=$(D_MPS)_$(Lx)x$(Ly).png",fig)








