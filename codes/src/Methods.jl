

function showdomain(M::AbstractTensorMap)
    @show codomain(M),domain(M)
end

function showQuantSweep(lsQ::Vector;name::String = "Quantity")
    for (iq,q) in enumerate(lsQ)
        println("$name\t$iq\t$q")
    end
end

function RandMPS(L::Int64,d::Int64,D::Int64)  
    
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
            iniMPS = U*S |> x -> x / norm(x)
            println("MPS initialized $(L-ii+1)/$L")
        end
        MPS[1] = permute(iniMPS,(),(1,2))
        println("MPS initialized $L/$L")
        println("MPS totally initialized")
        MPS
    end

    return MPS
end

function Quant1(ψ::Vector,Q::AbstractTensorMap,D_MPS::Int64)

    C = Vector{AbstractTensorMap}(undef,L)
    lsQi = Vector{Float64}(undef,L)
    C[1] = ψ[1]
    for i in 1:L

        if i ==1
            Qi = @tensor C[i][1,2]*Q[1,3]*C[i]'[3,2]
        elseif i == L 
            Qi = @tensor C[i][1,2]*Q[2,4]*C[i]'[1,4]
        else
            Qi = @tensor C[i][1,2,3]*Q[2,4]*C[i]'[1,4,3]
        end
        
        if imag(Qi) > 1e-5
            @error "non hermit"
        else
            lsQi[i] = real(Qi)
        end

        if i < L
            C[i:i+1] = RightMove(ψ[i+1],C[i],D_MPS)
        end
    end
    return lsQi
end
