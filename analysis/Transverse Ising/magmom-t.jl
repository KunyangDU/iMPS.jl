using CairoMakie,JLD2,TensorKit,LaTeXStrings

#include("../../calculations/iMPS/iMPS.jl")
#include("../../calculations/trans Ising/model.jl")

L = 9
J = -1.0
h = -0.8
D_MPS = 2^6

Mi = load("../codes/examples/Transverse Ising/data/tdvp/Impur_Mi_D=$(D_MPS)_L=$(L)_J=$(J)_h=$(h).jld2")["Mi"]
lst = load("../codes/examples/Transverse Ising/data/tdvp/Impur_lst_D=$(D_MPS)_L=$(L)_J=$(J)_h=$(h).jld2")["lst"]

width,height = 0.7 .* (600,200)

fig = Figure()
ax = Axis(fig[1,1],
xticks = 1:13,
xlabel = L"\text{site}\ i",
ylabel = L"t/J",
title = "TransIsing Dynamics\nL=$(L) J=$(J) h=$(h) D=$(D_MPS) (hz=0.01)",
titlealign = :left,
width = width,height = height)

hm = heatmap!(ax,1:L,lst,Mi',
colormap = :bwr,
#colorrange = (0,0.25)
)
#ylims!(ax,0,8)

Colorbar(fig[1,2],hm,
label = L"\langle σ^z \rangle")
resize_to_layout!(fig)
display(fig)

save("Transverse Ising/figures/Impur_Mi_D=$(D_MPS)_L=$(L)_J=$(J).pdf",fig)



