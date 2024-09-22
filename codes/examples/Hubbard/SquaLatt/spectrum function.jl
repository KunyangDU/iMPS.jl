using TensorKit,JLD2,FiniteLattices
include("../model.jl")
include("../../../src/iMPS.jl")


Lx = 6
Ly = 1
Latt = YCSqua(Lx,Ly)
μ= 0.0
t = 1

d = 4

D_MPS = 20


TruncErr = 1e-2
MaxIter = 20
ipath = [0 2*pi;0 0]
#ipath = [0 pi pi 0;0 0 pi 0]
kvecpath = vrange(ipath;eachstep = 2*size(Latt))
kr = pathlength(kvecpath)
ωscale = 10
lsω = collect(-1.2:0.2:1.2)*ωscale

for U in [0,8]


ψ = load("examples/Hubbard/data/$(Lx)x$(Ly)/ψ_D=$(D_MPS)_$(Lx)x$(Ly)_t=$(t)_U=$(U)_μ=$(μ).jld2")["ψ"]
Latt = load("examples/Hubbard/data/$(Lx)x$(Ly)/Latt_$(Lx)x$(Ly).jld2")["Latt"]
lsE = load("examples/Hubbard/data/$(Lx)x$(Ly)/lsE_D=$(D_MPS)_$(Lx)x$(Ly)_t=$(t)_U=$(U)_μ=$(μ).jld2")["lsE"]

H = HamMPO(Latt;μ=μ,U=U,t=t,d=d)
#= ψ = RandMPS(Lx*Ly)
ψ,lsE = sweepDMRG2(ψ,H,Nsweep,LanczosLevel,D_MPS) =#


Skω = zeros(length(kr),length(lsω),2)

for (ki,kv) in enumerate(kvecpath |> x -> collect.(eachcol(x)))
    println("--------------$(ki)----------------")

    for (oi,opr) in enumerate([aup,F*adown])

        ck = KOprMPO(Latt,opr,kv,-1;d=d,string=F)
        ckdagg = KOprMPO(Latt,collect(opr'),kv,1;d=d,string=F)

        Gk = GreenFuncRet(ψ,H,lsE[end],ck,ckdagg,lsω,D_MPS;TruncErr=TruncErr,MaxIter=MaxIter,τ=1e-2)
        Skω[ki,:,oi] = -imag.(Gk) / pi
    end

end

@save "examples/Hubbard/data/$(Lx)x$(Ly)/Skω_D=$(D_MPS)_$(Lx)x$(Ly)_t=$(t)_U=$(U).jld2" Skω
@save "examples/Hubbard/data/$(Lx)x$(Ly)/lsω_$(Lx)x$(Ly).jld2" lsω
@save "examples/Hubbard/data/$(Lx)x$(Ly)/kr_$(Lx)x$(Ly).jld2" kr
GC.gc()
end

