include("../../../../src/iMPS.jl")
include("../../model.jl")

Lx = 8
Ly = 4
d = 2
t=1

D_MPS = 2^3


Latt = load("examples/Free Fermion/data/$(Lx)x$(Ly)/Latt_$(Lx)x$(Ly).jld2")["Latt"]
lsμ = load("examples/Free Fermion/data/$(Lx)x$(Ly)/lsμ_$(Lx)x$(Ly).jld2")["lsμ"]
H = HamMPO(Latt)
ipath = [0 pi pi 0;0 0 pi 0]
#ipath = [-pi pi;0 0]

kvecpath = vrange(ipath;eachstep = size(Latt)-1)
kr = pathlength(kvecpath)

OccuN = zeros(length(kr),length(lsμ))
χk = zeros(length(kr),length(lsμ)-1)
totalNk = zeros(length(lsμ))

for (ki,kv) in enumerate(kvecpath |> x -> collect.(eachcol(x)))
    ck = CKMPO(Latt,kv)
    ckdagg = CKdaggMPO(Latt,kv)
    for (iμ,μ) in enumerate(lsμ)
        println("--------------$(ki),$(iμ)----------------")
        ψ = load("examples/Free Fermion/data/$(Lx)x$(Ly)/ψ_D=$(D_MPS)_$(Lx)x$(Ly)_t=$(t)_μ=$(μ).jld2")["ψ"]

        ckψ = VariContract(ck,ψ,D_MPS)
        #nkψ = VariContract(ckdagg,ckψ,D_MPS) ckdagg有问题
    
        OccuN[ki,iμ] = ApproxReal(InnerProd(ckψ,ckψ))
    end
end

for i in 1:length(lsμ)-1
    χk[:,i] = (OccuN[:,i + 1] - OccuN[:,i])/(lsμ[i + 1]-lsμ[i])/Lx/Ly
end

for i in 1:length(lsμ)
    totalNk[i] = (sum([OccuN[ki,i]*(kr[ki+1]-kr[ki]) for ki in eachindex(kr)[1:end-1]]) + sum([OccuN[ki,i]*(kr[ki]-kr[ki-1]) for ki in eachindex(kr)[2:end]]))/2
end

@save "examples/Free Fermion/data/$(Lx)x$(Ly)/χk_D=$(D_MPS)_$(Lx)x$(Ly).jld2" χk
@save "examples/Free Fermion/data/$(Lx)x$(Ly)/χkr_D=$(D_MPS)_$(Lx)x$(Ly).jld2" kr
@save "examples/Free Fermion/data/$(Lx)x$(Ly)/totalNk_D=$(D_MPS)_$(Lx)x$(Ly).jld2" totalNk





