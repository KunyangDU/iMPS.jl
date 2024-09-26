using CairoMakie,JLD2,TensorKit,LaTeXStrings

include("../../src/methods.jl")
include("../model.jl")

Lx = 21
Ly = 1
J = 1
h = 0.3
D_MPS = 2^5
t = 40*J
name = "Center"

lsObsDict = load("../codes/examples/Transverse Ising/data/lsObsDict_$(Lx)x$(Ly)_D_MPS=$(D_MPS)_J=$(J)_h=$(h)_time=$(t)_name=$(name).jld2")["lsObsDict"]
lst = load("../codes/examples/Transverse Ising/data/lst_$(Lx)x$(Ly)_D_MPS=$(D_MPS)_J=$(J)_h=$(h)_time=$(t)_name=$(name).jld2")["lst"]

Sz = ObsMat1(lsObsDict,"Sz")

width,height = 0.7 .* (600,200)
centralize
fig = Figure()
ax = Axis(fig[1,1],
xticks = 1:Lx,
xlabel = L"\text{site}\ i",
ylabel = L"tJ",
title = "TransIsing Dynamics\nL=$(Lx), h=$(h / J)J, D=$(D_MPS)",
titlealign = :left,
width = width,height = height)

hm = heatmap!(ax,1:Lx,lst,Sz,
colormap = :bwr,
#colorrange = (-0.3,0.3)
)
xlims!(ax,0.5,Lx + 0.5)
#ylims!(ax,0,8)

Colorbar(fig[1,2],hm,
label = L"\langle S^z \rangle")

vm = MaxGroupV(J,h)
tv = range(extrema(lst)...,10)
lines!(ax,vm*tv .+ (Lx-1)/2 .+ 1,tv,color = :white,linestyle = :dash)
lines!(ax,-vm*tv .+ (Lx-1)/2 .+ 1,tv,color = :white,linestyle = :dash)
text!(14, 8, text = L"x=x_0 + v^g_{max}t", align = (:left,:center),color = :white)

resize_to_layout!(fig)
display(fig)

save("Transverse Ising/figures/quench_$(Lx)x$(Ly)_D_MPS=$(D_MPS)_J=$(J)_h=$(h)_time=$(t).png",fig)



