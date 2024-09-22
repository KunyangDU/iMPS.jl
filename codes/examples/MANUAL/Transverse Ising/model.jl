
function HamMPO(L::Int64;J::Number=1,h::Number=0,hz::Number = 1e-2)

    I = diagm(ones(2))
    I0 = zeros(2,2)
    σx = [0 1;1 0]
    σz = [1 0;0 -1]

    d = 2
    D = 3

    MPO = Vector{AbstractTensorMap{ComplexSpace,2,2}}(undef, L)

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

function MagmomMPO(L::Int64)

    I = diagm(ones(2))
    I0 = zeros(2,2)
    σz = [1 0;0 -1]

    d = 2
    D = 2

    MPO = Vector{AbstractTensorMap{ComplexSpace,2,2}}(undef, L)

    idt = ℂ^1
    phys = (ℂ^d)'
    bond = ℂ^D
    for i in 1:L
        if i == 1
            data = reshape([σz I],d,1,d,D)
            M = BlockMPO(data,phys,idt,phys,bond)
        elseif i == L
            data = reshape([I;σz],d,D,d,1)
            M = BlockMPO(data,phys,bond,phys,idt)
        else
            data = reshape([
                I I0
                σz I
            ],d,D,d,D)
            M = BlockMPO(data,phys,bond,phys,bond)
        end
        MPO[i] = M
    end
    println("MPO constructed")

    return MPO
end

function LocalMagmomMPO(L::Int64,site::Int64;h::Number=0,t::Number=0)

    I = diagm(ones(2))
    I0 = zeros(2,2)
    σx = [0 1;1 0]
    σz = [1 0;0 -1]
    Σ = [I0 for _ in 1:L]
    Σ[site] = exp(-1im*h*σx*t)'*σz*exp(-1im*h*σx*t)

    d = 2
    D = 2

    MPO = Vector{AbstractTensorMap{ComplexSpace,2,2}}(undef, L)

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

function LocalMagmomMPO(d::Int64=2,data::Array=[1 0;0 -1])
    return Opr1(d,data)
end

function FerroMPS(L::Int64,state::String)
    # suppose center at leftmost

    if state == "FM"
        Σ = [reshape([1 0],1,2,1) for _ in 1:L]
    elseif state == "AFM"
        if iseven(L)
            Σ = repeat([reshape([1 0],1,2,1), reshape([0 1],1,2,1)], div(L,2))
        else
            Σ = vcat([reshape([1 0],1,2,1)],repeat([reshape([0 1],1,2,1), reshape([1 0],1,2,1)], div(L,2)))
        end
    else
        @error "state doesn't exist"
    end

    MPS = Vector{Union{AbstractTensorMap{ComplexSpace,1,2},AbstractTensorMap{ComplexSpace,0,3}}}(undef,L)

    bond = ℂ^1
    phys = (ℂ^2)'

    MPS[1] = Tensor(Σ[1],bond' ⊗ phys' ⊗ bond') |> x -> permute(x,(),(1,2,3)) / norm(x)
    
    for i in 2:L
        MPS[i] = let
            TensorMap(Σ[i],bond,phys ⊗ bond) |> x -> x / norm(x)
        end
    end
    
    return MPS
end

function ImpurMPS(L::Int64,site::Int64)
    # suppose center at leftmost

    Σ = [reshape([0 1],1,2,1) for _ in 1:L]
    Σ[site] = reshape([1 0],1,2,1)

    MPS = Vector{Union{AbstractTensorMap{ComplexSpace,1,2},AbstractTensorMap{ComplexSpace,0,3}}}(undef,L)

    bond = ℂ^1
    phys = (ℂ^d)'

    MPS[1] = Tensor(Σ[1],bond' ⊗ phys' ⊗ bond') |> x -> permute(x,(),(1,2,3)) / norm(x)
    
    for i in 2:L
        MPS[i] = let
            TensorMap(Σ[i],bond,phys ⊗ bond) |> x -> x / norm(x)
        end
    end
    
    return MPS
end

