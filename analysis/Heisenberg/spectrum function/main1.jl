using CairoMakie,JLD2,TensorKit,LaTeXStrings
include("../../src/MPSanalysis.jl")
include("../model.jl")

Lx = 8
Ly = 1

Jx = 1
Jy = Jx
Jz = Jx

d = 2
D_MPS = 2^3

Skω = load("../codes/examples/Heisenberg/data/$(Lx)x$(Ly)/Skω_D=$(D_MPS)_$(Lx)x$(Ly)_J=($((Jx,Jy,Jz))).jld2")["Skω"]
lsω = load("../codes/examples/Heisenberg/data/$(Lx)x$(Ly)/lsω_$(Lx)x$(Ly).jld2")["lsω"]
kr = load("../codes/examples/Heisenberg/data/$(Lx)x$(Ly)/kr_$(Lx)x$(Ly).jld2")["kr"]
lsE = load("../codes/examples/Heisenberg/data/$(Lx)x$(Ly)/lsE_D=$(D_MPS)_$(Lx)x$(Ly)_J=($((Jx,Jy,Jz))).jld2")["lsE"]

width,height = 1 .* (300,200)

xtickvalues = [0,pi/2,pi,3*pi/2,2*pi]
xticklabels = [L"0",L"\pi/2",L"\pi",L"3\pi/2",L"2\pi"]

fig = Figure()
ax = Axis(fig[1,1],
ylabel = L"\omega/t/4",
xlabel = L"k_r",
xticks = (xtickvalues,xticklabels),
title = "$(Lx)x$(Ly) Squa Heisenberg: J=$((Jx,Jy,Jz)), D=$(D_MPS)",
titlealign = :left,
yticks = 0:0.5:4,
width = width,height = height)
#hm = heatmap!(ax,kr,lsω / 4 / pi,Skω)
hm = heatmap!(ax,kr,lsω / 4 ,Skω)
Colorbar(fig[1,2],hm,
label = L"S(k,\ \omega)")


ipath = [0 2*pi;0 0]
kvecpath = vrange(ipath;eachstep = 100) 
theokr = pathlength(kvecpath)

# for antiferro, spectrum lies between πJ|sin k/2| and πJ|sin k|/2
#= upedge = @. abs(Jx*sin(theokr/2)) *pi
downedge = @. abs(Jx*sin(theokr)/2) *pi
lines!(ax,theokr,upedge,color = :red,linestyle = :dash)
lines!(ax,theokr,downedge,color = :red,linestyle = :dash) =#

# for ferro, spectrum lies at 
edge = @. 1*Jx*(1-cos(theokr))
lines!(ax,theokr,edge,color = :red,linestyle = :dash)
#ylims!(ax,0,3)

resize_to_layout!(fig)
display(fig)

save("Heisenberg/figures/Skω_D=$(D_MPS)_$(Lx)x$(Ly)_J=($((Jx,Jy,Jz))).pdf",fig)

