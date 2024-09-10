
function RandMPS(L::Int64;d::Int64=2)
    # suppose center at leftmost
    MPS = Vector{AbstractTensorMap}(undef,L)

    bond = ℂ^1
    phys = (ℂ^d)'

    MPS[1] = Tensor(rand(1,2,1),bond' ⊗ phys' ⊗ bond') |> x -> permute(x,(),(1,2,3)) / norm(x)
    
    for i in 2:L
        MPS[i] = let
            TensorMap(rand(1,2,1),bond,phys ⊗ bond) |> x -> x / norm(x)
        end
    end
    
    return MPS
end

function ApproxReal(Qi::Number;tol::Float64=1e-5)
    imag(Qi) <= tol && return real(Qi)
    @error "non hermit"
end


function LocalMerge(ψ1::AbstractTensorMap{ComplexSpace,1,2},
    ψ2::AbstractTensorMap{ComplexSpace,0,3})

    @tensor ψm[-1,-2,-3,-4] ≔ ψ1[1,-1,-2] * ψ2[1,-3,-4]
    return permute(ψm,(),(1,2,3,4))
end

function LocalMerge(ψ1::AbstractTensorMap{ComplexSpace,0,3},
    ψ2::AbstractTensorMap{ComplexSpace,1,2})

    @tensor ψm[-1,-2,-3,-4] ≔ ψ1[-1,-2,1] * ψ2[1,-3,-4]
    return permute(ψm,(),(1,2,3,4))
end



function InnerProd(ψ₁::Vector,ψ₂::Vector)
    innerprod = ψ₁[1]*ψ₂[1]'
    return innerprod[1]
end

