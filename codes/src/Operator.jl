

function BlockMPO(data::Array,
    cophys::ElementarySpace,cobond::ElementarySpace,
    phys::ElementarySpace,bond::ElementarySpace)
    
    return permute(TensorMap(data,cophys ⊗ cobond',phys ⊗ bond' ),(4,1),(2,3))
end

function BlockMPO(data::Array,
    cophys::ElementarySpace,phys::ElementarySpace)
    return TensorMap(data,cophys ,phys)
end

function EffHam(ψ::Vector,H::Vector,site::Int64)

    HR = RightEnv(ψ,H,site)
    HL = LeftEnv(ψ,H,site)

    return EffHam(H[site],HL,HR)
end

function EffHam(EnvL::AbstractTensorMap{ComplexSpace,2,1},EnvR::AbstractTensorMap{ComplexSpace,1,2})
    @tensor effH[-1,-2,-3,-4] ≔ EnvL[-1,1,-3]*EnvR[-2,1,-4]
    return permute(effH,(1,2),(3,4))
end

function EffHam(Hi::AbstractTensorMap{ComplexSpace,2,2},EnvL::AbstractTensorMap{ComplexSpace,2,1},EnvR::AbstractTensorMap{ComplexSpace,1,2})
    @tensor effH[-1,-2,-3,-4,-5,-6] ≔ EnvL[-1,1,-4]*Hi[2,-2,1,-5]*EnvR[-3,2,-6]
    return permute(effH,(1,2,3),(4,5,6))
end

function EffHam(Hi::Vector{AbstractTensorMap},EnvL::AbstractTensorMap{ComplexSpace,2,1},EnvR::AbstractTensorMap{ComplexSpace,1,2})
    @tensor effH[-1,-2,-3,-4,-5,-6,-7,-8] ≔ EnvL[-1,2,-5]*Hi[1][1,-2,2,-6]*Hi[2][3,-3,1,-7]*EnvR[-4,3,-8]
    return permute(effH,(1,2,3,4),(5,6,7,8))
end



function EvolveOpr(Ham::AbstractTensorMap,τ::Float64)
    return exp(1im * Ham * τ)
end


function Apply(ψi::AbstractTensorMap{ComplexSpace,0,2},Opr::AbstractTensorMap{ComplexSpace,2,2})
    @tensor ψia[-1,-2] ≔ ψi[1,2]*Opr[1,2,-1,-2]
    return permute(ψia,(),(1,2))
end

function Apply(ψi::AbstractTensorMap{ComplexSpace,0,3},Opr::AbstractTensorMap{ComplexSpace,3,3})
    @tensor ψia[-1,-2,-3] ≔ ψi[1,3,2]*Opr[1,3,2,-1,-2,-3]
    return permute(ψia,(),(1,2,3))
end

function Apply(ψi::AbstractTensorMap{ComplexSpace,0,4},Opr::AbstractTensorMap{ComplexSpace,4,4})
    @tensor ψia[-1,-2,-3,-4] ≔ ψi[1,3,4,2]*Opr[1,3,4,2,-1,-2,-3,-4]
    return permute(ψia,(),(1,2,3,4))
end


function HamMPO(L::Int64;J::Number=1,h::Number=0,hz::Number = 1e-2)

    I = diagm(ones(2))
    I0 = zeros(2,2)
    σx = [0 1;1 0]
    σz = [1 0;0 -1]

    d = 2
    D = 3

    MPO = Vector{AbstractTensorMap}(undef, L)

    idt = ℂ^1
    phys = (ℂ^d)'
    bond = ℂ^D
    for i in 1:L
        if i == 1
            data = reshape([h*σx+hz*σz J*σz I],d,D,d,1)
            M = BlockMPO(data,phys,idt,phys,bond)
        elseif i == L
            data = reshape([I;σz;h*σx+hz*σz],d,D,d,1)
            M = BlockMPO(data,phys,bond,phys,idt)
        else
            data = reshape([
                I I0 I0
                σz I0 I0
                h*σx+hz*σz J*σz I
            ],d,D,d,D)
            M = BlockMPO(data,phys,bond,phys,bond)
        end
        MPO[i] = M
    end
    println("MPO constructed")

    return MPO
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

    MPO = Vector{AbstractTensorMap}(undef, L)

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
