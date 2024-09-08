
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

function Quant1(ψ::Vector,Q::AbstractTensorMap{ComplexSpace,1,1},D_MPS::Int64)

    L = length(ψ)

    ψcp = Vector{AbstractTensorMap}(undef,L)
    lsQi = Vector{Float64}(undef,L)
    
    ψcp[1] = deepcopy(ψ[1])
    for i in 1:L
#=         effQ = EffHam(Q,InitialLeftEnv(),InitialRightEnv())
        Qi = @tensor ψcp[i][1,5,2]*effQ[1,5,2,3,6,4]*ψcp[i]'[3,6,4] =#

        Qi = @tensor ψcp[i][1,3,2]*Q[3,4]*ψcp[i]'[1,4,2]

        lsQi[i] = ApproxReal(Qi)

        if i < L
            ψcp[i:i+1] = RightMove(ψ[i+1],ψcp[i],D_MPS)
        end
        
    end

    return lsQi
end

function QuantUniv(ψ::Vector,Q::Vector)
    effQ = EffHam(Q[1],LeftEnv(ψ,Q,1),RightEnv(ψ,Q,1))
    q = @tensor ψ[1][1,5,2]*effQ[1,5,2,3,6,4]*ψ[1]'[3,6,4]
    return ApproxReal(q)
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


