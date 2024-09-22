using CairoMakie,JLD2,TensorKit,LaTeXStrings
include("../../src/MPSanalysis.jl")
include("../model.jl")

Lx = 4
Ly = 4

t = 1
μ = 0.0

D_MPS = 2^3

Skω = load("../codes/examples/Free Fermion/data/$(Lx)x$(Ly)/Skω_D=$(D_MPS)_$(Lx)x$(Ly).jld2")["Skω"]
lsω = load("../codes/examples/Free Fermion/data/$(Lx)x$(Ly)/lsω_$(Lx)x$(Ly).jld2")["lsω"]
kr = load("../codes/examples/Free Fermion/data/$(Lx)x$(Ly)/kr_$(Lx)x$(Ly).jld2")["kr"]
lsE = load("../codes/examples/Free Fermion/data/$(Lx)x$(Ly)/lsE_D=$(D_MPS)_$(Lx)x$(Ly)_t=$(t)_μ=$(μ).jld2")["lsE"]

width,height = 1 .* (300,200)
xtickvalues = [0,pi,2*pi,(2+sqrt(2))*pi]
xticklabels = [L"(0,0)",L"(\pi,0)",L"(\pi,\pi)",L"(0,0)"]
#xtickvalues = [0,pi]
#xticklabels = [L"(0,0)",L"(\pi,0)"]
krpath = range(-pi,pi,length(kr))

fig = Figure()
ax = Axis(fig[1,1],
ylabel = L"\omega/t",
xlabel = L"k_r",
xticks = (xtickvalues,xticklabels),
title = "$(Lx)x$(Ly) Squa Free Fermion: D=$(D_MPS)",
yticks = -8:2:8,
width = width,height = height)
hm = heatmap!(ax,kr,lsω,Skω)
Colorbar(fig[1,2],hm,
label = L"S(k,\ \omega)")

ipath = [0 pi pi 0;0 0 pi 0]
kvecpath = vrange(ipath;eachstep = 100) 
theokr = pathlength(kvecpath)
band = [SquaBand(kv) for kv in collect.(eachcol(kvecpath))]
lines!(ax,theokr,band,color = :red,linestyle = :dash)
for xv in xtickvalues
    lines!(ax,[xv,xv],collect(extrema(lsω)),color = :black)
end

#= Nμ = load("../codes/examples/Free Fermion/data/$(Lx)x$(Ly)/Nμ_D=$(2^5)_$(Lx)x$(Ly).jld2")["Nμ"]
lsμ = load("../codes/examples/Free Fermion/data/$(Lx)x$(Ly)/lsμ_$(Lx)x$(Ly).jld2")["lsμ"]
totalNk = load("../codes/examples/Free Fermion/data/$(Lx)x$(Ly)/totalNk_D=$(D_MPS)_$(Lx)x$(Ly).jld2")["totalNk"]
axN = Axis(fig[1,3],
width = width,height = height)
scatterlines!(axN,lsμ,Nμ) =#


resize_to_layout!(fig)
display(fig)

save("Free Fermion/figures/Skω_D=$(D_MPS)_$(Lx)x$(Ly).pdf",fig)


