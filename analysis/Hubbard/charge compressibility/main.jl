using CairoMakie,JLD2,TensorKit,LaTeXStrings
include("../../src/MPSanalysis.jl")
include("../model.jl")

Lx = 8
Ly = 1
U = 8

t = 1

D_MPS = 2^6

lsμ = load("../codes/examples/Hubbard/data/$(Lx)x$(Ly)/lsμ_$(Lx)x$(Ly)_U=$(U).jld2")["lsμ"]
Nμ = load("../codes/examples/Hubbard/data/$(Lx)x$(Ly)/Nμ_D=$(D_MPS)_$(Lx)x$(Ly)_U=$(U).jld2")["Nμ"]

nμ = Nμ / (Lx*Ly)
dμ = lsμ[2]-lsμ[1]

ind = 1:7:length(lsμ)
χ = diff(nμ[ind]) ./ diff(lsμ[ind])
centerμ = centralize(lsμ[ind])



width,height = 0.9 .* (300,200)

fig = Figure()

axμ = Axis(fig[1,1],
xlabel = L"n",
ylabel = L"μ",
title = "$(Lx)x$(Ly) Squa Hubbard: D=$(D_MPS)",
titlealign = :left,
width = width,height = height)
#ylims!(axμ,1.1.*extrema(lsμ)...)
scatterlines!(axμ,nμ[:,1],lsμ)


axχ = Axis(fig[1,2],
xlabel = L"χ",
xticks = 0:0.3:1,
width = 0.2*width,height = height)
hideydecorations!(axχ,grid = false)
scatterlines!(axχ,χ,centerμ)
xlims!(axχ,-0.1,0.6)

# U = 0, theoretical curve
#= ylims!(axχ,1.1.*extrema(lsμ)...)
theonμ = 0:0.01:1
theoμ = @. 2t*sin(pi*theonμ-pi/2)
theoχ = @. (1/pi) / sqrt((2*t)^2 - theoμ^2)
lines!(axμ,theonμ,theoμ,color = :red,linewidth = 2.0)
lines!(axχ,theoχ,theoμ,color = :red) =#

# U ≠ 0, band gap
lines!(axμ,[0,1],[3*U/4,3*U/4],color = :red,alpha=0.2)
lines!(axμ,[0,1],[U/4,U/4],color = :red,alpha=0.2)
text!(axμ,0.6, U/2, text = L"\Delta = U/2",
align = (:left,:center),color = :red)



resize_to_layout!(fig)

display(fig)

save("Hubbard/figures/χ_D=$(D_MPS)_$(Lx)x$(Ly)_U=$(U).pdf",fig)
