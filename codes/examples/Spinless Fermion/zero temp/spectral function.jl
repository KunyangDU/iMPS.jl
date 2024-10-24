include("../../../src/iMPS.jl")
include("../model.jl")

Lx = 16
Ly = 1
Latt = YCSqua(Lx,Ly)
t = 1
μ = 0.0

Nsweep = 5
LanczosLevel = 15
D_MPS = 2^5

H,D_MPO = compress(canonicalize(Hamiltonian(Latt;t=t,μ=μ)))

ψ = load("examples/Spinless Fermion/data/zero temp/ψ_$(Lx)x$(Ly)_D_MPS=$(D_MPS)_t=$(t)_μ=$(μ).jld2")["ψ"]
lsE = load("examples/Spinless Fermion/data/zero temp/lsE_$(Lx)x$(Ly)_D_MPS=$(D_MPS)_t=$(t)_μ=$(μ).jld2")["lsE"] 

τ = 5e-2
TruncErr = 1e-3
MaxIter = 50

ipath = [-pi pi;0 0]
#ipath = [0 pi pi 0;0 0 pi 0]
kvecpath = vrange(ipath;eachstep = size(Latt))
kr = pathlength(kvecpath)
ωscale = 5*pi # π for 1d
lsω = collect(-2*ωscale:0.05:2*ωscale)

Skω = zeros(length(kr),length(lsω))

for (ki,kv) in enumerate(kvecpath |> x -> collect.(eachcol(x)))
    println("--------------$(ki)----------------")

    kOprs,D_MPOs = map(x -> [mpo[x] for mpo in [ckMPO(Latt,kv),ckdagMPO(Latt,kv)]],1:2)
    Skω[ki,:] = SpecFuncFreq(ψ,kOprs...,H,lsE[end],D_MPS,LanczosLevel,lsω;τ=τ,TruncErr = TruncErr,MaxIter = MaxIter,ϵ=-1)

end

@save "examples/Spinless Fermion/data/zero temp/Skω_$(Lx)x$(Ly)_D_MPS=$(D_MPS)_t=$(t)_μ=$(μ).jld2" Skω
@save "examples/Spinless Fermion/data/zero temp/lsω_$(Lx)x$(Ly)_D_MPS=$(D_MPS)_t=$(t)_μ=$(μ).jld2" lsω
@save "examples/Spinless Fermion/data/zero temp/kr_$(Lx)x$(Ly)_D_MPS=$(D_MPS)_t=$(t)_μ=$(μ).jld2" kr



