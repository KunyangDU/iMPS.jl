

function HamMPO(Latt::AbstractLattice;
    t::Number=1,μ::Number=0,
    d::Int64=2)
    
    # JW transformation in Pauli series with convention: |0⟩ -> |↓⟩, |1⟩ -> |↑⟩
    # a⁺ -> σ⁺, a⁻ -> σ⁻
    # F = 1 - 2n = 1 - 2σ⁺σ⁻ -> -σz
    # https://zhuanlan.zhihu.com/p/386386413
    # https://itensor.org/docs.cgi?page=tutorials/fermions

    F = -[1 0;0 -1]
    a⁺ = [0 1;0 0]
    a = a⁺'

    maxd = FindMaxDist(neighbor(Latt))
    D_MPO = d*(2*maxd + 2)
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
                HM[k][(4*step-1)*d+1:(4*step+1)*d,(4*step-3)*d+1:(4*step-1)*d] = kron(diagm(ones(d)),F)
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

function HoppingMPO(HM::Vector,t::Number;
    d::Number = 2,D::Number = 7,
    F::Matrix = -[1 0;0 -1],a⁺::Matrix = [0 1;0 0],a::Matrix = a⁺')
    # suppose i < j
    Li = length(HM)
    hpHM = Vector{Matrix}(undef,Li)

    if Li == 2
        
    elseif Li == 3

    elseif Li == 4

    else
        @error "hopping doesn't exist"
    end


    return hpHM
end

function LocalNumMPO(L::Int64,site::Int64;)

    I = diagm(ones(2))
    I0 = zeros(2,2)
    a⁺ = [0 1;0 0]
    a = a⁺'
    n = a⁺*a

    Σ = [I0 for _ in 1:L]
    Σ[site] = n

    d = 2
    D = 2

    MPO = Vector{AbstractTensorMap}(undef, L)

    idt = ℂ^1
    phys = (ℂ^d)'
    bond = ℂ^D
    for i in 1:L
        if i == 1
            data = reshape([Σ[i] I],d,1,d,D)
            M = BlockMPO(data,phys,idt,phys,bond)
        elseif i == L
            data = reshape([I;Σ[i]],d,D,d,1)
            M = BlockMPO(data,phys,bond,phys,idt)
        else
            data = reshape([
                I I0
                Σ[i] I
            ],d,D,d,D)
            M = BlockMPO(data,phys,bond,phys,bond)
        end
        MPO[i] = M
    end
    println("MPO constructed")

    return MPO

end

function NMPO(L::Int64)

    I = diagm(ones(2))
    I0 = zeros(2,2)
    a⁺ = [0 1;0 0]
    a = a⁺'
    n = a⁺*a

    d = 2
    D = 2

    MPO = Vector{AbstractTensorMap}(undef, L)

    idt = ℂ^1
    phys = (ℂ^d)'
    bond = ℂ^D
    for i in 1:L
        if i == 1
            data = reshape([n I],d,D,d,1)
            M = BlockMPO(data,phys,idt,phys,bond)
        elseif i == L
            data = reshape([I;n],d,D,d,1)
            M = BlockMPO(data,phys,bond,phys,idt)
        else
            data = reshape([
                I I0
                n I
            ],d,D,d,D)
            M = BlockMPO(data,phys,bond,phys,bond)
        end
        MPO[i] = M
    end
    println("MPO constructed")

    return MPO
end


function LocalNumMPO(L::Int64,site::Int64;)

    I = diagm(ones(2))
    I0 = zeros(2,2)
    a⁺ = [0 1;0 0]
    a = a⁺'
    n = a⁺*a

    Σ = [I0 for _ in 1:L]
    Σ[site] = n

    d = 2
    D = 2

    MPO = Vector{AbstractTensorMap}(undef, L)

    idt = ℂ^1
    phys = (ℂ^d)'
    bond = ℂ^D
    for i in 1:L
        if i == 1
            data = reshape([Σ[i] I],d,1,d,D)
            M = BlockMPO(data,phys,idt,phys,bond)
        elseif i == L
            data = reshape([I;Σ[i]],d,D,d,1)
            M = BlockMPO(data,phys,bond,phys,idt)
        else
            data = reshape([
                I I0
                Σ[i] I
            ],d,D,d,D)
            M = BlockMPO(data,phys,bond,phys,bond)
        end
        MPO[i] = M
    end
    println("MPO constructed")

    return MPO

end

function LocalNMPO(d::Int64=2,data::Array=[1 0;0 0])
    return Opr1(d,data)
end



#= function FreeFermionDataSqua(t::Number = 1,μ::Number = 0)
    D = 8
    d = 2
    
    I = [1 0;0 1]
    F = -[1 0;0 -1]
    a⁺ = [0 1;0 0]
    a = a⁺'
    innerM = zeros(d*D,d*D)

    innerM[1:2,1:2] = I
    innerM[3:4,1:2] = a
    innerM[5:6,1:2] = a⁺
    innerM[15:16,1:2] = -μ*a⁺*a
    innerM[15:16,3:14] = repeat([t*a; t*a⁺],3)'
    innerM[15:16,15:16] = I
    innerM[7:14,3:10] = kron(diagm(ones(4)),F)

    return innerM
    
end =#

