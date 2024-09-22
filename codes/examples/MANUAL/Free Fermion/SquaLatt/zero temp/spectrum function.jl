include("../../../../src/iMPS.jl")
include("../../model.jl")



Lx = 6
Ly = 1
Latt = YCSqua(Lx,Ly)
d = 2
μ= 0.0
t=1

D_MPS = 2^3

ψ = load("examples/Free Fermion/data/$(Lx)x$(Ly)/ψ_D=$(D_MPS)_$(Lx)x$(Ly)_t=$(t)_μ=$(μ).jld2")["ψ"]
Latt = load("examples/Free Fermion/data/$(Lx)x$(Ly)/Latt_$(Lx)x$(Ly).jld2")["Latt"]
lsE = load("examples/Free Fermion/data/$(Lx)x$(Ly)/lsE_D=$(D_MPS)_$(Lx)x$(Ly)_t=$(t)_μ=$(μ).jld2")["lsE"]

H = HamMPO(Latt;μ=μ)
#= ψ = RandMPS(Lx*Ly)
ψ,lsE = sweepDMRG2(ψ,H,Nsweep,LanczosLevel,D_MPS) =#

TruncErr = 1e-2
MaxIter = 20
ipath = [-pi pi;0 0]
#ipath = [0 pi pi 0;0 0 pi 0]
kvecpath = vrange(ipath;eachstep = size(Latt))
kr = pathlength(kvecpath)
lsω = collect(-6.0:0.1:6.0)

Skω = Matrix{Float64}(undef,length(kr),length(lsω))

for (ki,kv) in enumerate(kvecpath |> x -> collect.(eachcol(x)))
    println("--------------$(ki)----------------")

    ck = KOprMPO(Latt,a,kv,-1;d=d,string = F)
    ckdagg = KOprMPO(Latt,a⁺,kv,1;d=d,string = F)

    Gk = GreenFuncRet(ψ,H,lsE[end],ck,ckdagg,lsω,D_MPS;TruncErr=TruncErr,MaxIter=MaxIter,τ=1e-2)
    Skω[ki,:] = -imag.(Gk) / pi
end

#= @save "examples/Free Fermion/data/$(Lx)x$(Ly)/Skω_D=$(D_MPS)_$(Lx)x$(Ly).jld2" Skω
@save "examples/Free Fermion/data/$(Lx)x$(Ly)/lsω_$(Lx)x$(Ly).jld2" lsω
@save "examples/Free Fermion/data/$(Lx)x$(Ly)/kr_$(Lx)x$(Ly).jld2" kr =#



