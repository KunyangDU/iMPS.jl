include("../../src/iMPS.jl")
include("model.jl")



Lx = 8
Ly = 1
J = -1
Latt = YCSqua(Lx,Ly)

LanczosLevel = 15
D_MPS = 2^4
TruncErr = 1e-3
MaxIter = 100


ψ = load("examples/Heisenberg/data/ψ_$(Lx)x$(Ly)_D_MPS=$(D_MPS)_J=$(J).jld2")["ψ"]
lsE = load("examples/Heisenberg/data/lsE_$(Lx)x$(Ly)_D_MPS=$(D_MPS)_J=$(J).jld2")["lsE"] 

H,D_MPO = compress(canonicalize(Hamiltonian(Latt;Jx=J,h=0)))

ipath = [0 2*pi;0 0]
#ipath = [0 pi pi 0;0 0 pi 0]
kvecpath = vrange(ipath;eachstep = size(Latt))
kr = pathlength(kvecpath)
ωscale = 5*pi # π for 1d
lsω = collect(0.0:0.05:2*ωscale)

Skω = zeros(length(kr),length(lsω))

for (ki,kv) in enumerate(kvecpath |> x -> collect.(eachcol(x)))
    println("--------------$(ki)----------------")
    SkOprs,D_MPOs = SkMPOs(Latt,kv)

    for SkOpr in SkOprs
        tempSqω = CorrFuncFreq(ψ,SkOpr,H,lsE[end],D_MPS,LanczosLevel,lsω;TruncErr=TruncErr,MaxIter=MaxIter,τ=1e-2,TR = true)
        Skω[ki,:] += ApproxReal.(tempSqω)
    end

end

@save "examples/Heisenberg/data/Skω_$(Lx)x$(Ly)_D_MPS=$(D_MPS)_J=$(J)_MaxIter=$(MaxIter).jld2" Skω
@save "examples/Heisenberg/data/lsω_$(Lx)x$(Ly)_D_MPS=$(D_MPS)_J=$(J)_MaxIter=$(MaxIter).jld2" lsω
@save "examples/Heisenberg/data/kr_$(Lx)x$(Ly)_D_MPS=$(D_MPS)_J=$(J)_MaxIter=$(MaxIter).jld2" kr

# 关联函数多个负号
# 算自由费米子的谱函数
# 下一步加对称性？


