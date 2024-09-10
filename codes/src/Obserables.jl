
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


function GreenFuncRet(ψ::AbstractTensorMap,H::Vector,E0::Float64,
    CK::Vector,CKdagg::Vector,lsE::Vector,D_MPS::Int64;
    τ::Number = 1/E0/5,TruncErr::Number = 1e-5,MaxIter::Int64=30)
    # calculate the G^{ret}(k,ω) with given cₖ,cₖ⁺
    # how to choose evolve step τ?
    # spectrum function is S(k,ω) = - Im G^{ret}(k,ω) / π

    Gk = Vector{ComplexF64}(undef,length(lsE))

    ψ₊E = Apply(ψ,CK)
    ψ₋E = Apply(ψ,CKdagg)

    lsGt₊,lst₊ = GreenFuncTDVP2(ψ₊E,H,τ,TruncErr,MaxIter,D_MPS)
    lsGt₋,lst₋ = GreenFuncTDVP2(ψ₋E,H,-τ,TruncErr,MaxIter,D_MPS)

    Gk = @. lsGt₋*exp(-1im*(lsE-E0)*lst₋) + lsGt₊*exp(-1im*(lsE-E0)*lst₊)
    
    return Gk
end




