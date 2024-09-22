using TensorKit,JLD2,FiniteLattices

include("../../../src/iMPS.jl")
include("../model.jl")

Lx = 4
Ly = 4

#= Jx = -1
Jy = Jx
Jz = Jx =#

d = 2
D_MPS = 2^3

TruncErr = 1e-2
MaxIter = 100

for Jx in [-1,1]
Jy = Jx
Jz = Jx

Latt = load("examples/Heisenberg/data/$(Lx)x$(Ly)/Latt_$(Lx)x$(Ly).jld2")["Latt"]
lsE = load("examples/Heisenberg/data/$(Lx)x$(Ly)/lsE_D=$(D_MPS)_$(Lx)x$(Ly)_J=($((Jx,Jy,Jz))).jld2")["lsE"]
ψ = load("examples/Heisenberg/data/$(Lx)x$(Ly)/ψ_D=$(D_MPS)_$(Lx)x$(Ly)_J=($((Jx,Jy,Jz))).jld2")["ψ"]
H = load("examples/Heisenberg/data/$(Lx)x$(Ly)/H_$(Lx)x$(Ly)_J=($((Jx,Jy,Jz))).jld2")["H"]
#= ψ = RandMPS(Lx*Ly)
ψ,lsE = sweepDMRG2(ψ,H,Nsweep,LanczosLevel,D_MPS) =#

#ipath = [0 2*pi;0 0]
ipath = [0 pi pi 0;0 0 pi 0]
kvecpath = vrange(ipath;eachstep = div(size(Latt),2))
kr = pathlength(kvecpath)
ωscale = 4*20 # 4*π for 1d, 4*8 for 2d
lsω = collect(0.0:0.01:1.2)*ωscale

Skω = zeros(length(kr),length(lsω))

σx = [0 1;1 0]
σy = [0 -1im;1im 0]
σz = [1 0;0 -1]
σs = [σx, σy, σz]

for (ki,kv) in enumerate(kvecpath |> x -> collect.(eachcol(x)))
    println("--------------$(ki)----------------")

    for σi in σs
        σk = KOprMPO(Latt,σi,kv,-1)
    
        totalSqω = DynStrucFac(ψ,H,lsE[end],σk,lsω,D_MPS;TruncErr=TruncErr,MaxIter=MaxIter,τ=1e-2,)
        Skω[ki,:] += totalSqω
    end

    #@show kv,Skω[ki,:]
    
end

@save "examples/Heisenberg/data/$(Lx)x$(Ly)/Skω_D=$(D_MPS)_$(Lx)x$(Ly)_J=($((Jx,Jy,Jz))).jld2" Skω
@save "examples/Heisenberg/data/$(Lx)x$(Ly)/lsω_$(Lx)x$(Ly).jld2" lsω
@save "examples/Heisenberg/data/$(Lx)x$(Ly)/kr_$(Lx)x$(Ly).jld2" kr
end