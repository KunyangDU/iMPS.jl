function RightMove(nextψ::AbstractTensorMap,eigenAi::AbstractTensorMap,D_MPS::Int64)
    Σ,V = RightSVD(eigenAi,D_MPS)
    return collect(RightMerge(Σ,V,nextψ))
end

function LeftMove(nextψ::AbstractTensorMap,eigenAi::AbstractTensorMap,D_MPS::Int64)
    Σ,V = LeftSVD(eigenAi,D_MPS)
    return collect(LeftMerge(Σ,V,nextψ))
end

# SVD

function RightSVD(eigenAi::AbstractTensorMap{ComplexSpace, 0, 2},D_MPS::Int64)
    U,S,V = tsvd(eigenAi,(2,),(1,);trunc = truncdim(D_MPS))
    return permute(U*S,(),(1,2)),V
end

function RightSVD(eigenAi::AbstractTensorMap{ComplexSpace, 0, 3},D_MPS::Int64)
    U,S,V = tsvd(eigenAi,(3,),(1,2);trunc = truncdim(D_MPS))
    return permute(U*S,(),(1,2)),V
end

function LeftSVD(eigenAi::AbstractTensorMap{ComplexSpace,0,2},D_MPS::Int64)
    U,S,V = tsvd(eigenAi,(1,),(2,);trunc = truncdim(D_MPS))
    return permute(U*S,(),(1,2)),V
end

function LeftSVD(eigenAi::AbstractTensorMap{ComplexSpace,0,3},D_MPS::Int64)
    U,S,V = tsvd(eigenAi,(1,),(2,3);trunc = truncdim(D_MPS))
    return permute(U*S,(),(1,2)),V
end

# MERGE

function RightMerge(Σ::AbstractTensorMap{ComplexSpace, 0, 2},V::AbstractTensorMap{ComplexSpace, 1, 1},nextψ::AbstractTensorMap{ComplexSpace, 1, 2})
    @tensor tempMPS[-1,-2,-3] ≔ Σ[1,-1]*nextψ[1,-2,-3]
    return V,permute(tempMPS,(),(1,2,3))
end

function RightMerge(Σ::AbstractTensorMap{ComplexSpace, 0, 2},V::AbstractTensorMap{ComplexSpace, 1, 2},nextψ::AbstractTensorMap{ComplexSpace, 1, 1})
    @tensor tempMPS[-1,-2] ≔ Σ[-1,1]*nextψ[1,-2]
    return V,permute(tempMPS,(),(1,2))
end

function RightMerge(Σ::AbstractTensorMap{ComplexSpace, 0, 2},V::AbstractTensorMap{ComplexSpace, 1, 2},nextψ::AbstractTensorMap{ComplexSpace, 1, 2})
    @tensor tempMPS[-1,-2,-3] ≔ Σ[-1,1]*nextψ[1,-2,-3]
    return V,permute(tempMPS,(),(1,2,3))
end

function LeftMerge(Σ::AbstractTensorMap{ComplexSpace, 0, 2},V::AbstractTensorMap{ComplexSpace, 1, 1},nextψ::AbstractTensorMap{ComplexSpace, 1, 2})
    @tensor tempMPS[-1,-2,-3] ≔ Σ[1,-1]*nextψ[1,-2,-3]
    return permute(tempMPS,(),(1,2,3)),V
end

function LeftMerge(Σ::AbstractTensorMap{ComplexSpace, 0, 2},V::AbstractTensorMap{ComplexSpace, 1, 2},nextψ::AbstractTensorMap{ComplexSpace, 1, 1})
    @tensor tempMPS[-1,-2] ≔ Σ[1,-1]*nextψ[1,-2]
    return permute(tempMPS,(),(2,1)),V
end

function LeftMerge(Σ::AbstractTensorMap{ComplexSpace, 0, 2},V::AbstractTensorMap{ComplexSpace, 1, 2},nextψ::AbstractTensorMap{ComplexSpace, 1, 2})
    @tensor tempMPS[-1,-2,-3] ≔ Σ[1,-1]*nextψ[1,-2,-3]
    return permute(tempMPS,(),(1,2,3)),V
end

# more stable
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


