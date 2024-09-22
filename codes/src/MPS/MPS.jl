
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


