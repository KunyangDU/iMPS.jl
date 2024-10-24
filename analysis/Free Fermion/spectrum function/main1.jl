using CairoMakie,JLD2,TensorKit,LaTeXStrings,FiniteLattices
include("../../src/MPSanalysis.jl")
include("../model.jl")

Lx = 8
Ly = 1

t = 1
μ = 0.0

D_MPS = 2^3

Skω = load("../codes/examples/Spinless Fermion/data/zero temp/Skω_$(Lx)x$(Ly)_D_MPS=$(D_MPS)_t=$(t)_μ=$(μ).jld2")["Skω"]
lsω = load("../codes/examples/Spinless Fermion/data/zero temp/lsω_$(Lx)x$(Ly)_D_MPS=$(D_MPS)_t=$(t)_μ=$(μ).jld2")["lsω"]
kr = load("../codes/examples/Spinless Fermion/data/zero temp/kr_$(Lx)x$(Ly)_D_MPS=$(D_MPS)_t=$(t)_μ=$(μ).jld2")["kr"]
lsE = load("../codes/examples/Spinless Fermion/data/zero temp/lsE_$(Lx)x$(Ly)_D_MPS=$(D_MPS)_t=$(t)_μ=$(μ).jld2")["lsE"]

width,height = 1 .* (300,200)

xtickvalues = [0,pi/2,pi,3*pi/2,2*pi]
xticklabels = [L"-\pi",L"-\pi/2",L"0",L"\pi/2",L"\pi"]

fig = Figure()
ax = Axis(fig[1,1],
ylabel = L"\omega/t",
xlabel = L"k_r",
xticks = (xtickvalues,xticklabels),
title = "$(Lx)x$(Ly) Squa Free Fermion: D=$(D_MPS)",
yticks = -4:2:4,
width = width,height = height)
hm = heatmap!(ax,kr,lsω,Skω)
Colorbar(fig[1,2],hm,
label = L"S(k,\ \omega)")

ipath = [-pi pi;0 0]
kvecpath = vrange(ipath;eachstep = 100) 
theokr = pathlength(kvecpath)
band = [SquaBand(kv)-2 for kv in collect.(eachcol(kvecpath))]
lines!(ax,theokr,band,color = :red,linestyle = :dash)
for xv in xtickvalues
    lines!(ax,[xv,xv],collect(extrema(lsω)),color = :black)
end

ylims!(ax,-4,4)

resize_to_layout!(fig)
display(fig)

save("Free Fermion/figures/Skω_D=$(D_MPS)_$(Lx)x$(Ly).pdf",fig)


