using CairoMakie,JLD2,TensorKit,LaTeXStrings
include("../../src/MPSanalysis.jl")
include("../model.jl")

Lx = 8
Ly = 1

t = 1
μ = 0.0
U = 8

D_MPS = 2^5

Skω = load("../codes/examples/AUTO/Hubbard/data/Skω_$(Lx)x$(Ly)_D_MPS=$(D_MPS)_t=$(t)_μ=$(μ)_U=$(U).jld2")["Skω"]
lsω = load("../codes/examples/AUTO/Hubbard/data/lsω_$(Lx)x$(Ly)_D_MPS=$(D_MPS)_t=$(t)_μ=$(μ)_U=$(U).jld2")["lsω"]
kr = load("../codes/examples/AUTO/Hubbard/data/kr_$(Lx)x$(Ly)_D_MPS=$(D_MPS)_t=$(t)_μ=$(μ)_U=$(U).jld2")["kr"]
lsE = load("../codes/examples/AUTO/Hubbard/data/lsE_$(Lx)x$(Ly)_D_MPS=$(D_MPS)_t=$(t)_μ=$(μ)_U=$(U).jld2")["lsE"]

width,height = 1 .* (300,200)

xtickvalues = [0,pi/2,pi,3*pi/2,2*pi]
xticklabels = [L"-\pi",L"-\pi/2",L"0",L"\pi/2",L"\pi"]

fig = Figure()
ax = Axis(fig[1,1],
ylabel = L"\omega/t",
xlabel = L"k_r",
xticks = (xtickvalues,xticklabels),
title = "$(Lx)x$(Ly) Squa Hubbard: D=$(D_MPS)",
yticks = -8:2:8,
width = width,height = height)
hm = heatmap!(ax,kr,lsω,Skω[:,end:-1:1,2])
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

ylims!(ax,-12,12)

resize_to_layout!(fig)
display(fig)

save("Hubbard/figures/Skω_$(Lx)x$(Ly)_D_MPS=$(D_MPS)_t=$(t)_μ=$(μ)_U=$(U).jld2",fig)


