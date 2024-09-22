

function BlockMPO(data::Array,
    cophys::ElementarySpace,cobond::ElementarySpace,
    phys::ElementarySpace,bond::ElementarySpace)
    
    return permute(TensorMap(data,cophys ⊗ cobond',phys ⊗ bond' ),(4,1),(2,3))
end

function BlockMPO(data::Array,
    cophys::ElementarySpace,phys::ElementarySpace)
    return TensorMap(data,cophys ,phys)
end


function Opr1(d::Int64,data::Array)
    phys = (ℂ^d)'
    return BlockMPO(data,phys,phys)
end

function UnivMPO(L::Int64,Opr::Matrix)
    Σ = [Opr for _ in 1:L]
    return LocalMPO(L,Σ)
end

function LocalMPO(L::Int64,Opr::Matrix{T},site::Int64) where T
    d = size(Opr)[1]
    Σ = [zeros(T,d,d) for _ in 1:L]
    Σ[site] = Opr
    return LocalMPO(L,Σ)
end

function LocalMPO(L::Int64,lsOpr::Vector)

    d = size(lsOpr[1])[1]

    D = 2
    I = diagm(ones(d))
    I0 = zeros(d,d)

    MPO = Vector{AbstractTensorMap{ComplexSpace,2,2}}(undef, L)

    idt = ℂ^1
    phys = (ℂ^d)'
    bond = ℂ^D
    for i in 1:L
        if i == 1
            data = reshape([lsOpr[i] I],d,1,d,D)
            M = BlockMPO(data,phys,idt,phys,bond)
        elseif i == L
            data = reshape([I;lsOpr[i]],d,D,d,1)
            M = BlockMPO(data,phys,bond,phys,idt)
        else
            data = reshape([
                I I0
                lsOpr[i] I
            ],d,D,d,D)
            M = BlockMPO(data,phys,bond,phys,bond)
        end
        MPO[i] = M
    end
    println("MPO constructed")

    return MPO
    
end

#= function Apply(MPS::Vector,MPO::Vector)
    

    return finalMPS
end =#

function KOprMPO(Latt::AbstractLattice,Opr::Matrix,kv::Vector,sign::Int64;
    d::Int64 = 2,D::Int64 = 2,
    string::Matrix = diagm(ones(size(Opr)[1])))

    # take h.c. of Opr automatically
    # sign denotes the +- sign in the exp(...)

    L = size(Latt)

    I = diagm(ones(d))
    I0 = zeros(d,d)

    MPO = Vector{AbstractTensorMap}(undef, L)

    idt = ℂ^1
    phys = (ℂ^d)'
    bond = ℂ^D
    for i in 1:L
        if i == 1
            data = reshape([Opr'*exp(sign*1im*dot(kv,coordinate(Latt,i)))/sqrt(L) string],d,1,d,D)
            M = BlockMPO(data,phys,idt,phys,bond)
        elseif i == L
            data = reshape([I;Opr'*exp(sign*1im*dot(kv,coordinate(Latt,i)))/sqrt(L)],d,D,d,1)
            M = BlockMPO(data,phys,bond,phys,idt)
        else
            data = reshape([
                I I0
                Opr'*exp(sign*1im*dot(kv,coordinate(Latt,i)))/sqrt(L) string
            ],d,D,d,D)
            M = BlockMPO(data,phys,bond,phys,bond)
        end
        MPO[i] = M
    end
    println("MPO constructed")

    return MPO
    
end

function RandMPO(L::Int64,d::Int64)

    MPO = Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}}(undef,L)

    idt = ℂ^1
    phys = (ℂ^d)'


    MPO[1] = TensorMap(rand(ComplexF64,1,d,1,d),phys,idt ⊗ phys ⊗ idt)
    
    for i in 2:L
        MPO[i] = let
            TensorMap(rand(ComplexF64,1,d,1,d),phys ⊗ idt,phys ⊗ idt)
        end
    end
    println("MPO constructed")

    return MPO
    
end


function EvolveOpr(Ham::AbstractTensorMap,τ::Float64)
    return exp(1im * Ham * τ)
end


function IdentityMPO(L::Int64,d::Int64)
    
    MPO = Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}}(undef,L)

    idt = ℂ^1
    phys = (ℂ^d)'

    MPO[1] = TensorMap(diagm(ones(d)),phys,idt ⊗ phys ⊗ idt)
    
    for i in 2:L
        MPO[i] = let
            TensorMap(diagm(ones(d)),phys ⊗ idt,phys ⊗ idt)
        end
    end
    println("MPO constructed")

    return MPO
end

