using TensorKit,JLD2,FiniteLattices,CairoMakie

include("model.jl")
include("../../src/iMPS.jl")

#= function diagm(pairs...)
    L = length(pairs[1][2]) + pairs[1][1]
    mat = zeros(L,L)
    for pair in pairs
        mat += diagm(pair)
    end
    return mat
end =#

Lx = 6
Ly = 6
Latt = YCSqua(Lx,Ly)

ipath = [0 pi pi 0;0 0 pi 0]
kvecpath = vrange(ipath;eachstep = 2*size(Latt)-1)
kr = pathlength(kvecpath)

kvecpathg = kdivide(kvecpath,size(ipath)[2]-1)
krg = kdivide(kr,size(ipath)[2]-1)
