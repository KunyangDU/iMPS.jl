using CairoMakie,JLD2,TensorKit,LaTeXStrings

#include("../../calculations/iMPS/iMPS.jl")
#include("../../calculations/trans Ising/model.jl")

L = 13
J = -1.0
h = 0.5
D_MPS = 2^4

#DiffMi = load("trans Ising/data/tdvp/DiffMi_D=$(D_MPS)_L=$(L)_J=$(J)_h=$(h).jld2")["DiffMi"]
Mi = load("trans Ising/data/tdvp/Mi_D=$(D_MPS)_L=$(L)_J=$(J)_h=$(h).jld2")["DiffMi"]
lst = load("trans Ising/data/tdvp/lst_D=$(D_MPS)_L=$(L)_J=$(J)_h=$(h).jld2")["lst"]

width,height = 0.9 .* (300,200)

fig = Figure()
ax = Axis(fig[1,1],
xticks = 1:13,
xlabel = L"\text{site}\ i",
ylabel = L"t",
title = "TransIsing Dynamics\nL=$(L) J=$(J) D=$(D_MPS) (hz=0.01)",
titlealign = :left,
width = width,height = height)

hm = heatmap!(ax,1:L,lst,Mi',
colormap = :bwr,
#colorrange = (0,0.25)
)
#ylims!(ax,0,3)

Colorbar(fig[1,2],hm,
label = L"\langle σ^z \rangle")
resize_to_layout!(fig)
display(fig)

save("trans Ising/figures/DiffMi_D=$(D_MPS)_L=$(L)_J=$(J).pdf",fig)
