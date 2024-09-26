using CairoMakie,JLD2,TensorKit,LaTeXStrings,FiniteLattices
include("../../src/MPSanalysis.jl")
include("../model.jl")

Lx = 16
Ly = 1
U = 0

Latt = YCSqua(Lx,Ly)

t = 1

D_MPS = 2^6

lsμ = load("../codes/examples/AUTO/Hubbard/data/lsμ_$(Lx)x$(Ly)_D_MPS=$(D_MPS)_t=$(t)_U=$(U).jld2")["lsμ"]
Nμ = zeros(length(lsμ),2)
for (μi,μ) in enumerate(lsμ)
    ObsDict = load("../codes/examples/AUTO/Hubbard/data/ObsDict_$(Lx)x$(Ly)_D_MPS=$(D_MPS)_t=$(t)_μ=$(μ)_U=$(U).jld2")["ObsDict"]
    Nμ[μi,1],Nμ[μi,2] = map(x -> sum([ObsDict[x][(i,)] for i in 1:size(Latt)]),("nup","ndown"))
end

nμ = Nμ / (Lx*Ly)
dμ = lsμ[2]-lsμ[1]

ind = 1:1:length(lsμ)
χ = diff(nμ[ind,1]) ./ diff(lsμ[ind]) + diff(nμ[ind,2]) ./ diff(lsμ[ind])
centerμ = centralize(lsμ[ind])



width,height = 0.9 .* (300,200)

fig = Figure()

axμ = Axis(fig[1,1],
xlabel = L"n",
ylabel = L"μ",
title = "$(Lx)x$(Ly) Squa Hubbard: D=$(D_MPS),U=$(U)",
titlealign = :left,
width = width,height = height)
#ylims!(axμ,1.1.*extrema(lsμ)...)
scatterlines!(axμ,nμ[:,1] + nμ[:,2],lsμ)


axχ = Axis(fig[1,2],
xlabel = L"χ",
xticks = 0:0.2:1,
width = 0.2*width,height = height)
hideydecorations!(axχ,grid = false)
scatterlines!(axχ,χ,centerμ)
#xlims!(axχ,-0.1,1.2)

# U = 0, theoretical curve
ylims!(axχ,1.1.*extrema(lsμ)...)
theonμ = 0:0.01:2
theoμ = @. 2t*sin(pi*theonμ/2-pi/2)
theoχ = @. (2/pi) / sqrt((2*t)^2 - theoμ^2)
lines!(axμ,theonμ,theoμ,color = :red,linewidth = 2.0)
lines!(axχ,theoχ,theoμ,color = :red)

# U ≠ 0, band gap
#= lines!(axμ,[0,2],[3*U/4,3*U/4],color = :red,alpha=0.2)
lines!(axμ,[0,2],[U/4,U/4],color = :red,alpha=0.2)
text!(axμ,1.2, U/2, text = L"\Delta = U/2",
align = (:left,:center),color = :red) =#



resize_to_layout!(fig)

display(fig)

save("Hubbard/figures/χ__$(Lx)x$(Ly)_D_MPS=$(D_MPS)_t=$(t)_U=$(U).pdf",fig)
