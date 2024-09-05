
function initialMPS(L::Int64,d::Int64,D::Int64)  
    
    MPS = let 
        iniMPS = Tensor(rand, ComplexF64, ⊗([ℂ^d for i in 1:L]...)) |> x -> x / norm(x)
        MPS = Vector{AbstractTensorMap}(undef, L)
        for ii in L:-1:2
            if ii == L
                U,S,V = tsvd(iniMPS,Tuple.((1:ii-1,ii:L))...;trunc = truncdim(D))
                MPS[ii] = V
            else
                U,S,V = tsvd(iniMPS,Tuple.((1:ii-1,ii:ii+1))...;trunc = truncdim(D))
                MPS[ii] = V
            end
            iniMPS = U*S
            println("MPS initialized $(L-ii+1)/$L")
        end
        MPS[1] = permute(iniMPS,(),(1,2))
        println("MPS initialized $L/$L")
        println("MPS totally initialized")
        MPS
    end

    return MPS
end


