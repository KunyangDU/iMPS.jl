
function RandMPS(L::Int64;d::Int64=2)
    # suppose center at leftmost
    MPS = Vector{Union{AbstractTensorMap{ComplexSpace,0,3},AbstractTensorMap{ComplexSpace,1,2}}}(undef,L)

    bond = ℂ^1
    phys = (ℂ^d)'

    MPS[1] = Tensor(rand(ComplexF64,1,d,1),bond' ⊗ phys' ⊗ bond') |> x -> permute(x,(),(1,2,3)) / norm(x)
    
    for i in 2:L
        MPS[i] = let
            TensorMap(rand(1,d,1),bond,phys ⊗ bond) |> x -> x / norm(x)
        end
    end
    
    return MPS
end


function ManualMPS(Latt::AbstractLattice,Σ::Vector)

    L = size(Latt)

    MPS = Vector{Union{AbstractTensorMap{ComplexSpace,1,2},AbstractTensorMap{ComplexSpace,0,3}}}(undef,L)

    d = length(Σ[1])
    
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

function ManualMPS(Latt::AbstractLattice,name::String;d::Int64=2)

    L = size(Latt)

    if name == "AFM"
        if iseven(L)
            Σ = repeat([[1 0], [0 1]], div(L,2))
        else
            Σ = vcat([[1 0]],repeat([[0 1], [1 0]], div(L,2)))
        end
    elseif name == "FM"
        Σ = [[1 0] for _ in 1:L]
    elseif name == "Center"
        Σ = [[0 1] for _ in 1:L]
        Σ[round(Int64,(L-1)/2) + 1] = [1 0]
    else
        @error "Name not defined!"
    end

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


