include("../../src/iMPS.jl")
include("model.jl")


Lx = 8
Ly = 1
Latt = YCSqua(Lx,Ly)
t = 1
U = 8
d = 4
μ = 0.0

lsμ = (U/4 - 4):0.5:(3*U/4 + 4)

D_MPS = 2^5
LanczosLevel = 15
Nsweep = 5

H,D_MPO = compress(canonicalize(Hamiltonian(Latt;t=t,μ=μ,U=U)))

ψ = load("examples/Hubbard/data/ψ_$(Lx)x$(Ly)_D_MPS=$(D_MPS)_t=$(t)_μ=$(μ)_U=$(U).jld2")["ψ"]
lsE = load("examples/Hubbard/data/lsE_$(Lx)x$(Ly)_D_MPS=$(D_MPS)_t=$(t)_μ=$(μ)_U=$(U).jld2")["lsE"] 
τ = 5e-2
TruncErr = 1e-2
MaxIter = 100

ipath = [-pi pi;0 0]
#ipath = [0 pi pi 0;0 0 pi 0]
kvecpath = vrange(ipath;eachstep = 2*size(Latt))
kr = pathlength(kvecpath)
ωscale = 5*pi # π for 1d
lsω = collect(-2*ωscale:0.05:2*ωscale)

Skω = zeros(length(kr),length(lsω),2)

for (ki,kv) in enumerate(kvecpath |> x -> collect.(eachcol(x)))
    println("--------------$(ki)----------------")
    spinKOprs,D_MPOs = cMPO(Latt,kv)
    for (si,kOprs) in enumerate(spinKOprs)
        Skω[ki,:,si] = SpecFuncFreq(ψ,kOprs...,H,lsE[end],D_MPS,LanczosLevel,lsω;τ=τ,TruncErr = TruncErr,MaxIter = MaxIter,ϵ=-1)
    end
end

@save "examples/Hubbard/data/Skω_$(Lx)x$(Ly)_D_MPS=$(D_MPS)_t=$(t)_μ=$(μ)_U=$(U).jld2" Skω
@save "examples/Hubbard/data/lsω_$(Lx)x$(Ly)_D_MPS=$(D_MPS)_t=$(t)_μ=$(μ)_U=$(U).jld2" lsω
@save "examples/Hubbard/data/kr_$(Lx)x$(Ly)_D_MPS=$(D_MPS)_t=$(t)_μ=$(μ)_U=$(U).jld2" kr

