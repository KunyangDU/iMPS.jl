

function HamMPO(L::Int64;Jx::Number=1,Jy::Number=Jx,Jz::Number=Jx,hz::Number = 1e-2)

    I = diagm(ones(2))
    I0 = zeros(2,2)
    σx = [0 1;1 0]
    σy = [0 -1im;1im 0]
    σz = [1 0;0 -1]

    d = 2
    D = 5

    MPO = Vector{AbstractTensorMap}(undef, L)

    idt = ℂ^1
    phys = (ℂ^d)'
    bond = ℂ^D
    for i in 1:L
        if i == 1
            data = reshape([-hz*σz -Jx*σx -Jy*σy -Jz*σz I],d,D,d,1)
            M = BlockMPO(data,phys,idt,phys,bond)
        elseif i == L
            data = reshape([I;σx;σy;σz;-hz*σz],d,D,d,1)
            M = BlockMPO(data,phys,bond,phys,idt)
        else
            data = reshape([
                I I0 I0 I0 I0
                σx I0 I0 I0 I0
                σy I0 I0 I0 I0
                σz I0 I0 I0 I0
                -hz*σz -Jx*σx -Jy*σy -Jz*σz I
            ],d,D,d,D)
            M = BlockMPO(data,phys,bond,phys,bond)
        end
        MPO[i] = M
    end
    println("MPO constructed")

    return MPO
end



function HamMPO(Latt::AbstractLattice;
    Jx::Number=1,Jy::Number=Jx,Jz::Number=Jx,hz::Number = 1e-2,
    d::Int64=2)
    
    # https://zhuanlan.zhihu.com/p/391850098

    I = diagm(ones(2))
    I0 = zeros(2,2)
    σx = [0 1;1 0]
    σy = [0 -1im;1im 0]
    σz = [1 0;0 -1]

    d = 2

    maxd = FindMaxDist(neighbor(Latt))
    D_MPO = 3*maxd + 2
    HM = InitHamMatri(size(Latt),d,D_MPO)

    # onsite pin field
    mode = zeros(D_MPO,D_MPO)
    mode[end,1] = 1
    HM[1] = HM[1] .+ kron(mode[end:end,:],-hz*σz)
    for i in eachindex(HM)[2:end-1]
        HM[i] = HM[i] .+ kron(mode,-hz*σz)
    end
    HM[end] = HM[end] .+ kron(mode[:,1:1],-hz*σz)

    # NN hopping

    pairs = neighbor(Latt)
    # 如果hopping参数不一样，那么每个参数对应的hopping应该多一个态
    #for iLatt in size(Latt):-1:2
    for iLatt in size(Latt):-1:2
        hps = FindPair(pairs,2,iLatt)
        for (i,j) in hps
            dist = j-i
            HM[j][d .+ (1:d),1:d] = σx
            HM[j][2*d .+ (1:d),1:d] = σy'
            HM[j][3*d .+ (1:d),1:d] = σz

            HM[i][end-d+1:end,(3*dist-2)*d .+ (1:d)] = -Jx*σx
            HM[i][end-d+1:end,(3*dist-1)*d .+ (1:d)] = -Jy*σy'
            HM[i][end-d+1:end,3*dist*d .+ (1:d)] = -Jz*σz

            for k in j-1:-1:i+1
                step = j-k
                #HM[k][(4*step-1)*d+1:(4*step+1)*d,(4*step-3)*d+1:(4*step-1)*d] = kron(diagm(ones(d)),F)
                HM[k][(3*step+1)*d .+ (1:3*d),(3*step-2)*d .+ (1:3*d)] = kron(diagm(ones(3)),I)
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
    HM = Vector{Matrix{ComplexF64}}(undef,L)

    mode = kron(diagm(vcat([1],zeros(D_MPO-2),[1])),diagm(ones(d)))

    HM[1] = mode[end-d+1:end,:]
    HM[end] = mode[:,1:d]
    for i in eachindex(HM)[2:end-1]
        HM[i] = deepcopy(mode)
    end

    return HM
end





