
function RightSVD(nextψ::AbstractTensorMap,Ev::AbstractTensorMap,site::Int64,D_MPS::Int64)

    if site == 1
        U,S,V = tsvd(Ev,(2,),(1,);trunc = truncdim(D_MPS))
        @tensor tempMPS[-1,-2,-3] ≔ permute(U*S,(),(1,2))[1,-1]*nextψ[1,-2,-3]
        nextMPS = permute(tempMPS,(),(1,2,3))
    elseif site == length(ψ)-1
        U,S,V = tsvd(Ev,(3,),(1,2);trunc = truncdim(D_MPS))
        @tensor tempMPS[-1,-2] ≔ permute(U*S,(),(1,2))[-1,1]*nextψ[1,-2]
        nextMPS = permute(tempMPS,(),(1,2))
    else
        U,S,V = tsvd(Ev,(3,),(1,2);trunc = truncdim(D_MPS))
        @tensor tempMPS[-1,-2,-3] ≔ permute(U*S,(),(1,2))[-1,1]*nextψ[1,-2,-3]
        nextMPS = permute(tempMPS,(),(1,2,3))
    end
    thisMPS = V

    return [thisMPS,nextMPS]

end

function LeftSVD(nextψ::AbstractTensorMap,Ev::AbstractTensorMap,site::Int64,D_MPS::Int64)

    if site == length(ψ)
        U,S,V = tsvd(Ev,(1,),(2,);trunc = truncdim(D_MPS))
        @tensor tempMPS[-1,-2,-3] ≔ permute(U*S,(),(1,2))[1,-1]*nextψ[1,-2,-3]
        nextMPS = permute(tempMPS,(),(1,2,3))
    elseif site == 2
        U,S,V = tsvd(Ev,(1,),(2,3);trunc = truncdim(D_MPS))
        @tensor tempMPS[-1,-2] ≔ permute(U*S,(),(1,2))[1,-1]*nextψ[1,-2]
        nextMPS = permute(tempMPS,(),(2,1))
    else
        U,S,V = tsvd(Ev,(1,),(2,3);trunc = truncdim(D_MPS))
        @tensor tempMPS[-1,-2,-3] ≔ permute(U*S,(),(1,2))[1,-1]*nextψ[1,-2,-3]
        nextMPS = permute(tempMPS,(),(1,2,3))
    end
    thisMPS = V

    return [nextMPS,thisMPS]

end

function DirectedSVD(ψ::Vector,eigenTM::AbstractTensorMap,site::Int64,direction::String,D_MPS::Int64)
    if direction == "right"

        if site == 1
            U,S,V = tsvd(eigenTM,(2,),(1,);trunc = truncdim(D_MPS))
            @tensor tempMPS[-1,-2,-3] ≔ permute(U*S,(),(1,2))[1,-1]*ψ[site+1][1,-2,-3]
            nextMPS = permute(tempMPS,(),(1,2,3))
        elseif site == length(ψ)-1
            U,S,V = tsvd(eigenTM,(3,),(1,2);trunc = truncdim(D_MPS))
            @tensor tempMPS[-1,-2] ≔ permute(U*S,(),(1,2))[-1,1]*ψ[site+1][1,-2]
            nextMPS = permute(tempMPS,(),(1,2))
        else
            U,S,V = tsvd(eigenTM,(3,),(1,2);trunc = truncdim(D_MPS))
            @tensor tempMPS[-1,-2,-3] ≔ permute(U*S,(),(1,2))[-1,1]*ψ[site+1][1,-2,-3]
            nextMPS = permute(tempMPS,(),(1,2,3))
        end
        thisMPS = V

        MPSs = [thisMPS,nextMPS]
    elseif direction == "left"

        if site == length(ψ)
            U,S,V = tsvd(eigenTM,(1,),(2,);trunc = truncdim(D_MPS))
            @tensor tempMPS[-1,-2,-3] ≔ permute(U*S,(),(1,2))[1,-1]*ψ[site-1][1,-2,-3]
            nextMPS = permute(tempMPS,(),(1,2,3))
        elseif site == 2
            U,S,V = tsvd(eigenTM,(1,),(2,3);trunc = truncdim(D_MPS))
            @tensor tempMPS[-1,-2] ≔ permute(U*S,(),(1,2))[1,-1]*ψ[site-1][1,-2]
            nextMPS = permute(tempMPS,(),(2,1))
        else
            U,S,V = tsvd(eigenTM,(1,),(2,3);trunc = truncdim(D_MPS))
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
