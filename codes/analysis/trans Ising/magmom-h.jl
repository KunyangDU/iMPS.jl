using CairoMakie,JLD2,TensorKit,LaTeXStrings

#include("../../calculations/iMPS/iMPS.jl")
#include("../../calculations/trans Ising/model.jl")

L = 12
J = -1.0
D_MPS = 2^4

@load "trans Ising/data/dmrg/lsh_D=$(D_MPS)_L=$(L)_J=$(J).jld2" lsh
@load "trans Ising/data/dmrg/lsMz_D=$(D_MPS)_L=$(L)_J=$(J).jld2" lsMz
@load "trans Ising/data/dmrg/lsEg_D=$(D_MPS)_L=$(L)_J=$(J).jld2" lsEg


width,height = 0.9 .* (450,110)


fig = Figure()
axMz = Axis(fig[1,1],
ylabel = L"Mz",
title = "Transverse Ising: L=$(L) J=$(J) D=$(D_MPS) (hz=0.01)",
titlealign = :left,
width = width,height = height)
ylims!(-13,1)

scatterlines!(axMz,lsh,lsMz)

hidexdecorations!(axMz,grid = false)


axEg = Axis(fig[2,1],
xlabel = L"h",
ylabel = L"Eg",
width = width,height = height)
scatterlines!(axEg,lsh,lsEg)

resize_to_layout!(fig)

display(fig)

save("trans Ising/figures/D=$(D_MPS)_L=$(L)_J=$(J).pdf",fig)
