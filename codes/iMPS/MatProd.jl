
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

function RightEnv(ψ::Vector,H::Vector,site::Int64=1)
    @tensor HR[-1,-2,-3] ≔ ψ[L][-1,1]*H[L][1,2,-2]*ψ[L]'[2,-3]

    HR = permute(HR,(1,),(2,3))
    for iL in L-1:-1:site+1
        @tensor B[-1,-2,-3] ≔ ψ[iL][-1,1,4]*H[iL][1,5,2,-2]*ψ[iL]'[2,6,-3]*HR[4,5,6]
        HR = permute(B,(1,),(2,3))
    end

    return HR
end

function LeftEnv(ψ::Vector,H::Vector,site::Int64)

    @tensor HL[-1,-2,-3] ≔ ψ[1][-1,1]*H[1][1,-2,2]*ψ[1]'[2,-3]
    HL = permute(HL,(1,2),(3,))
    for iL in 2:site-1
        @tensor B[-1,-2,-3] ≔ ψ[iL][-1,4,1]*H[iL][1,-2,2,5]*ψ[iL]'[6,2,-3]*HL[4,5,6]
        HL = permute(B,(1,2),(3,))
    end

    return HL
end

function reduHam(ψ::Vector,H::Vector,site::Int64)

    # 缩并顺序可以优化

    if site == 1
        HR = RightEnv(ψ,H,site)
        @tensor reduH[-1,-2,-3,-4] ≔ H[site][-1,1,-3]*HR[-2,1,-4]
        reduH = permute(reduH,(1,2),(3,4))
    elseif site == length(ψ)
        HL = LeftEnv(ψ,H,site)
        @tensor reduH[-1,-2,-3,-4] ≔ HL[-1,1,-3]*H[site][-2,-4,1]
        reduH = permute(reduH,(1,2),(3,4))
    else
        HR = RightEnv(ψ,H,site)
        HL = LeftEnv(ψ,H,site)
        @tensor reduH[-1,-2,-3,-4,-5,-6] ≔ HL[-1,1,-4]*H[site][-2,2,-5,1]*HR[-3,2,-6]
        reduH = permute(reduH,(1,2,3),(4,5,6))
    end

    return reduH
end


function orientSVD(ψ::Vector,eigenTM::AbstractTensorMap,site::Int64,direction::String)
    if direction == "right"

        if site == 1
            U,S,V = tsvd(eigenTM,(2,),(1,))
            @tensor tempMPS[-1,-2,-3] ≔ permute(U*S,(),(1,2))[1,-1]*ψ[site+1][1,-2,-3]
            nextMPS = permute(tempMPS,(),(1,2,3))
        elseif site == length(ψ)-1
            U,S,V = tsvd(eigenTM,(3,),(1,2))
            @tensor tempMPS[-1,-2] ≔ permute(U*S,(),(1,2))[-1,1]*ψ[site+1][1,-2]
            nextMPS = permute(tempMPS,(),(1,2))
        else
            U,S,V = tsvd(eigenTM,(3,),(1,2))
            @tensor tempMPS[-1,-2,-3] ≔ permute(U*S,(),(1,2))[-1,1]*ψ[site+1][1,-2,-3]
            nextMPS = permute(tempMPS,(),(1,2,3))
        end
        thisMPS = V

        MPSs = [thisMPS,nextMPS]
    elseif direction == "left"

        if site == length(ψ)
            U,S,V = tsvd(eigenTM,(1,),(2,))
            @tensor tempMPS[-1,-2,-3] ≔ permute(U*S,(),(1,2))[1,-1]*ψ[site-1][1,-2,-3]
            nextMPS = permute(tempMPS,(),(1,2,3))
        elseif site == 2
            U,S,V = tsvd(eigenTM,(1,),(2,3))
            @tensor tempMPS[-1,-2] ≔ permute(U*S,(),(1,2))[1,-1]*ψ[site-1][1,-2]
            nextMPS = permute(tempMPS,(),(1,2))
        else
            U,S,V = tsvd(eigenTM,(1,),(2,3))
            @tensor tempMPS[-1,-2,-3] ≔ permute(U*S,(),(1,2))[1,-1]*ψ[site-1][1,-2,-3]
            nextMPS = permute(tempMPS,(),(1,2,3))
        end
        thisMPS = V

        MPSs = [nextMPS,thisMPS]
    else
        @error "key word error"
    end

    return MPSs
end
