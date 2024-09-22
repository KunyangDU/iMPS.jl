using CairoMakie,JLD2,TensorKit,LaTeXStrings

#include("../../calculations/iMPS/iMPS.jl")
#include("../../calculations/trans Ising/model.jl")

L = 12
J = -1.0
h = 0.01
D_MPS = 2^5

@load "../codes/examples/Transverse Ising/data/dmrg/lsh_D=$(D_MPS)_L=$(L)_J=$(J).jld2" lsh
@load "../codes/examples/Transverse Ising/data/dmrg/lsMz_D=$(D_MPS)_L=$(L)_J=$(J).jld2" lsMz
@load "../codes/examples/Transverse Ising/data/dmrg/lsEg_D=$(D_MPS)_L=$(L)_J=$(J).jld2" lsEg


width,height = 0.9 .* (450,110)


fig = Figure()
axMz = Axis(fig[1,1],
ylabel = L"\langle \sigma^z \rangle",
title = "Transverse Ising: L=$(L) J=$(J) hz =$(h) D=$(D_MPS) (hz=0.01)",
titlealign = :left,
width = width,height = height)
ylims!(-1.2,0.2)

scatterlines!(axMz,lsh,lsMz / L)


theomz = @. -(1-min(1,abs(lsh/J)^2))^(1/8)
lines!(axMz,lsh,theomz,color = :red)


hidexdecorations!(axMz,grid = false)


axEg = Axis(fig[2,1],
xlabel = L"h",
ylabel = L"Eg",
width = width,height = height)
scatterlines!(axEg,lsh,lsEg)

resize_to_layout!(fig)

display(fig)

save("Transverse Ising/figures/D=$(D_MPS)_L=$(L)_J=$(J).pdf",fig)
