using CairoMakie


include("../model.jl")

J = 1.0
h = 1.0
L = 12

lsω = -1.0:0.2:5.0

width,height = 1 .* (300,200)

xtickvalues = [0,pi/2,pi,3*pi/2,2*pi]
xticklabels = [L"0",L"\pi/2",L"\pi",L"3\pi/2",L"2\pi"]

fig = Figure()
ax = Axis(fig[1,1],
ylabel = L"\omega/t",
xlabel = L"k",
xticks = (xtickvalues,xticklabels),
title = "L=$(L) TransIsing: D=$(D_MPS)",
yticks = -4:2:4,
width = width,height = height)

ipath = [0 2*pi;0 0]
kvecpath = vrange(ipath;eachstep = 100) 
theokr = pathlength(kvecpath)
band = [IsingBand(kv[1],J,h) for kv in collect.(eachcol(kvecpath))]
lines!(ax,theokr,band,color = :red,linestyle = :dash)
for xv in xtickvalues
    lines!(ax,[xv,xv],collect(extrema(lsω)),color = :black)
end

resize_to_layout!(fig)
display(fig)
