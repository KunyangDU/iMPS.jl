
function HamMPO(Latt::AbstractLattice;
    t::Number=1,μ::Number=0,
    d::Int64=2)
    
    # JW transformation in Pauli series with convention: |0⟩ -> |↓⟩, |1⟩ -> |↑⟩
    # cup = aup, cdown = F*adown
    # cup⁺ = aup⁺, cdown⁺ = adown⁺*F
    # https://zhuanlan.zhihu.com/p/386386413
    # https://itensor.org/docs.cgi?page=tutorials/fermions

    F = diagm([1,-1,-1,1])
    aup = diagm(-2 => [1,1])
    adown = diagm(2 => [1,1])

    maxd = FindMaxDist(neighbor(Latt))
    #D_MPO = d*(2*maxd + 2)
    D_MPO = 2*(maxd+1)
    HM = InitHamMatri(size(Latt),d,D_MPO)

    # onsite μ
    mode = zeros(D_MPO,D_MPO)
    mode[end,1] = 1
    HM[1] = HM[1] .+ kron(mode[end:end,:],-μ*a⁺*a)
    for i in eachindex(HM)[2:end-1]
        HM[i] = HM[i] .+ kron(mode,-μ*a⁺*a)
    end
    HM[end] = HM[end] .+ kron(mode[:,1:1],-μ*a⁺*a)

    # NN hopping

    pairs = neighbor(Latt)
    # 如果hopping参数不一样，那么每个参数对应的hopping应该多一个态
    #for iLatt in size(Latt):-1:2
    for iLatt in size(Latt):-1:2
        hps = FindPair(pairs,2,iLatt)
        for (i,j) in hps
            dist = j-i
            HM[j][d .+ (1:d),1:d] = a
            HM[j][2*d .+ (1:d),1:d] = a⁺

            HM[i][end-d+1:end,(2*dist-1)*d .+ (1:d)] = t*a⁺
            HM[i][end-d+1:end,2*dist*d .+ (1:d)] = t*a

            for k in j-1:-1:i+1
                step = j-k
                #HM[k][(4*step-1)*d+1:(4*step+1)*d,(4*step-3)*d+1:(4*step-1)*d] = kron(diagm(ones(d)),F)
                HM[k][(2*step+1)*d .+ (1:2*d),(2*step-1)*d .+ (1:2*d)] = kron(diagm(ones(d)),F)
            end
        end
    end

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



