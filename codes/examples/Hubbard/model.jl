global F = diagm([1,-1,-1,1])
global aup = diagm(-2 => [1,1])
global adown = diagm(-1 => [1,0,1])



function HamMPO(Latt::AbstractLattice;
    t::Number=1,U::Number=0,μ::Number=0,
    d::Int64=4)

    # JW transformation in Pauli series with convention: |0⟩ -> |↓⟩, |1⟩ -> |↑⟩
    # cup = aup, cdown = F*adown
    # cup⁺ = aup⁺, cdown⁺ = adown⁺*F
    # https://zhuanlan.zhihu.com/p/386386413
    # https://itensor.org/docs.cgi?page=tutorials/fermions

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
                HM[k][(4*step+1)*d .+ (1:4*d),(4*step-3)*d .+ (1:4*d)] = kron(diagm(ones(4)),F)
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

    MPO = Vector{AbstractTensorMap}(undef,size(Latt))
    
    idt = ℂ^1
    phys = (ℂ^d)'
    bond = ℂ^D_MPO
    MPO[1] = BlockMPO(reshape(HM[1],d,1,d,D_MPO),phys,idt,phys,bond)
    for i in eachindex(MPO)[2:end-1]
        MPO[i] = BlockMPO(reshape(HM[i],d,D_MPO,d,D_MPO),phys,bond,phys,bond)
    end
    MPO[end] = BlockMPO(reshape(HM[end],d,D_MPO,d,1),phys,bond,phys,idt)

    println("MPO constructed")

    return MPO
end


function InitHamMatri(L::Int64,d::Int64,D_MPO::Int64)
    HM = Vector{Matrix}(undef,L)

    mode = kron(diagm(vcat([1],zeros(D_MPO-2),[1])),diagm(ones(d)))

    HM[1] = mode[end-d+1:end,:]
    HM[end] = mode[:,1:d]
    for i in eachindex(HM)[2:end-1]
        HM[i] = deepcopy(mode)
    end

    return HM
end
