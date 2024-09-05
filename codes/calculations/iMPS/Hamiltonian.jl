

function reduHam(ψ::Vector,H::Vector,site::Int64)

    # 缩并顺序可以优化

    if site == 1
        HR = RightEnv(ψ,H,site)
        reduH = LeftMostReduHam(H[site],HR)
    elseif site == length(ψ)
        HL = LeftEnv(ψ,H,site)
        reduH = RightMostReduHam(H[site],HL)
    else
        HR = RightEnv(ψ,H,site)
        HL = LeftEnv(ψ,H,site)
        reduH = InnerReduHam(H[site],HL,HR)
    end

    return reduH
end

function LeftMostReduHam(Hi::AbstractTensorMap,HR::AbstractTensorMap)
    @tensor reduH[-1,-2,-3,-4] ≔ Hi[-1,1,-3]*HR[-2,1,-4]
    reduH = permute(reduH,(1,2),(3,4))
    return reduH
end

function RightMostReduHam(Hi::AbstractTensorMap,HL::AbstractTensorMap)
    @tensor reduH[-1,-2,-3,-4] ≔ HL[-1,1,-3]*Hi[-2,-4,1]
    reduH = permute(reduH,(1,2),(3,4))
    return reduH
end

function InnerReduHam(Hi::AbstractTensorMap,HL::AbstractTensorMap,HR::AbstractTensorMap)
    @tensor reduH[-1,-2,-3,-4,-5,-6] ≔ HL[-1,1,-4]*Hi[-2,2,-5,1]*HR[-3,2,-6]
    reduH = permute(reduH,(1,2,3),(4,5,6))
    return reduH
end

function RightSVD(ψ::Vector,Ev::AbstractTensorMap,site::Int64,D_MPS::Int64)

    if site == 1
        U,S,V = tsvd(Ev,(2,),(1,);trunc = truncdim(D_MPS))
        @tensor tempMPS[-1,-2,-3] ≔ permute(U*S,(),(1,2))[1,-1]*ψ[site+1][1,-2,-3]
        nextMPS = permute(tempMPS,(),(1,2,3))
    elseif site == length(ψ)-1
        U,S,V = tsvd(Ev,(3,),(1,2);trunc = truncdim(D_MPS))
        @tensor tempMPS[-1,-2] ≔ permute(U*S,(),(1,2))[-1,1]*ψ[site+1][1,-2]
        nextMPS = permute(tempMPS,(),(1,2))
    else
        U,S,V = tsvd(Ev,(3,),(1,2);trunc = truncdim(D_MPS))
        @tensor tempMPS[-1,-2,-3] ≔ permute(U*S,(),(1,2))[-1,1]*ψ[site+1][1,-2,-3]
        nextMPS = permute(tempMPS,(),(1,2,3))
    end
    thisMPS = V

    return [thisMPS,nextMPS]

end

function LeftSVD(ψ::Vector,Ev::AbstractTensorMap,site::Int64,D_MPS::Int64)

    if site == length(ψ)
        U,S,V = tsvd(Ev,(1,),(2,);trunc = truncdim(D_MPS))
        @tensor tempMPS[-1,-2,-3] ≔ permute(U*S,(),(1,2))[1,-1]*ψ[site-1][1,-2,-3]
        nextMPS = permute(tempMPS,(),(1,2,3))
    elseif site == 2
        U,S,V = tsvd(Ev,(1,),(2,3);trunc = truncdim(D_MPS))
        @tensor tempMPS[-1,-2] ≔ permute(U*S,(),(1,2))[1,-1]*ψ[site-1][1,-2]
        nextMPS = permute(tempMPS,(),(2,1))
    else
        U,S,V = tsvd(Ev,(1,),(2,3);trunc = truncdim(D_MPS))
        @tensor tempMPS[-1,-2,-3] ≔ permute(U*S,(),(1,2))[1,-1]*ψ[site-1][1,-2,-3]
        nextMPS = permute(tempMPS,(),(1,2,3))
    end
    thisMPS = V

    return [nextMPS,thisMPS]

end


function orientSVD(ψ::Vector,eigenTM::AbstractTensorMap,site::Int64,direction::String,D_MPS::Int64)
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
