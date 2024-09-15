using TensorKit,JLD2,FiniteLattices
include("../model.jl")
include("../../../src/iMPS.jl")

Lx = 4
Ly = 3
Latt = YCSqua(Lx,Ly)
@save "examples/Hubbard/data/$(Lx)x$(Ly)/Latt_$(Lx)x$(Ly).jld2" Latt

d = 4

lsμ = -4.0:0.2:4.0
@save "examples/Hubbard/data/$(Lx)x$(Ly)/lsμ_$(Lx)x$(Ly).jld2" lsμ

t = 1
U = 0
μ = 0.0

D_MPS = 2^3
maxd = FindMaxDist(neighbor(Latt))
D_MPO = d*(4*maxd + 2)

LanczosLevel = D_MPO*d
Nsweep = 3


H = HamMPO(Latt;t=t,μ=μ,U=U,d=d)

showBlockMPO(H[1])
1

F = diagm([1,-1,-1,1])
aup = diagm(-2 => [1,1])
adown = diagm(-1 => [1,0,1])

maxd = FindMaxDist(neighbor(Latt))
D_MPO = 4*maxd + 2
HM = InitHamMatri(size(Latt),d,D_MPO)

# onsite μ
mode = zeros(D_MPO,D_MPO)
mode[end,1] = 1
HM[1] = HM[1] .+ kron(mode[end:end,:],-μ*(aup'*F*aup + adown'*F*adown))
for i in eachindex(HM)[2:end-1]
    HM[i] = HM[i] .+ kron(mode,-μ*(aup'*F*aup + adown'*F*adown))
end
HM[end] = HM[end] .+ kron(mode[:,1:1],-μ*(aup'*F*aup + adown'*F*adown))

# NN hopping

pairs = neighbor(Latt)
# 如果hopping参数不一样，那么每个参数对应的hopping应该多一个态
#for iLatt in size(Latt):-1:2
# up down updagg downdagg
for iLatt in size(Latt):-1:2
    hps = FindPair(pairs,2,iLatt)
    for (i,j) in hps
        dist = j-i
        HM[j][d .+ (1:d),1:d]   = aup
        HM[j][2*d .+ (1:d),1:d] = F*adown
        HM[j][3*d .+ (1:d),1:d] = aup'
        HM[j][4*d .+ (1:d),1:d] = F*adown'

        HM[i][end-d+1:end,(4*dist-3)*d .+ (1:d)] =  t*aup'*F
        HM[i][end-d+1:end,(4*dist-2)*d .+ (1:d)] =  t*adown'
        HM[i][end-d+1:end,(4*dist-1)*d .+ (1:d)] = -t*aup*F
        HM[i][end-d+1:end,4*dist*d .+ (1:d)]     = -t*adown

        for k in j-1:-1:i+1
            step = j-k
            #HM[k][(4*step-1)*d+1:(4*step+1)*d,(4*step-3)*d+1:(4*step-1)*d] = kron(diagm(ones(d)),F)
            HM[k][(4*step+1)*d .+ (1:4*d),(4*step-3)*d .+ (1:4*d)] = diagm(ones(4*d))
        end
    end
end

# onsite U
mode = zeros(D_MPO,D_MPO)
mode[end,1] = 1
HM[1] = HM[1] .+ kron(mode[end:end,:],U*aup'*F*aup*adown'*F*adown)
for i in eachindex(HM)[2:end-1]
    HM[i] = HM[i] .+ kron(mode,U*aup'*F*aup*adown'*F*adown)
end
HM[end] = HM[end] .+ kron(mode[:,1:1],U*aup'*F*aup*adown'*F*adown)

HM[3][end-d+1:end,9*d .+ (1:4*d)]
