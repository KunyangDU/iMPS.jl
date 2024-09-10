function RightMove(nextψ::AbstractTensorMap{ComplexSpace, 1, 2},thisψ::AbstractTensorMap{ComplexSpace, 0, 3},D_MPS::Int64)
    V,Σ,truncerr = RightSVD(thisψ,D_MPS)
    return collect(RightMerge(Σ,V,nextψ)),truncerr
end

function LeftMove(nextψ::AbstractTensorMap{ComplexSpace, 1, 2},thisψ::AbstractTensorMap{ComplexSpace, 0, 3},D_MPS::Int64)
    Σ,V,truncerr = LeftSVD(thisψ,D_MPS)
    return collect(LeftMerge(Σ,V,nextψ)),truncerr
end

# SVD

function RightSVD(ψi::AbstractTensorMap{ComplexSpace, 0, 3},D_MPS::Int64)
    U,S,V,truncerr = tsvd(ψi,(3,),(1,2);trunc = truncdim(D_MPS))
    return V,permute(U*S,(),(1,2)),truncerr
end

function RightSVD(ψi::AbstractTensorMap{ComplexSpace, 0, 4},D_MPS::Int64)
    U,S,V,truncerr = tsvd(ψi,(3,4),(1,2);trunc = truncdim(D_MPS))
    return V,permute(U*S,(),(3,1,2)),truncerr
end

function LeftSVD(ψi::AbstractTensorMap{ComplexSpace,0,3},D_MPS::Int64)
    U,S,V,truncerr = tsvd(ψi,(1,),(2,3);trunc = truncdim(D_MPS))
    return permute(U*S,(),(1,2)),permute(V,(1,),(2,3)),truncerr
end

function LeftSVD(ψi::AbstractTensorMap{ComplexSpace,0,4},D_MPS::Int64)
    U,S,V,truncerr = tsvd(ψi,(1,2),(3,4);trunc = truncdim(D_MPS))
    return permute(U*S,(),(1,2,3)),permute(V,(1,),(2,3)),truncerr
end

# MERGE

function RightMerge(Σ::AbstractTensorMap{ComplexSpace, 0, 2},V::AbstractTensorMap{ComplexSpace, 1, 2},nextψ::AbstractTensorMap{ComplexSpace, 1, 2})
    @tensor tempMPS[-1,-2,-3] ≔ Σ[1,-1]*nextψ[1,-2,-3]
    return V,permute(tempMPS,(),(1,2,3))
end

function LeftMerge(Σ::AbstractTensorMap{ComplexSpace, 0, 2},V::AbstractTensorMap{ComplexSpace, 1, 2},nextψ::AbstractTensorMap{ComplexSpace, 1, 2})
    @tensor tempMPS[-1,-2,-3] ≔ Σ[1,-3]*nextψ[1,-1,-2]
    return permute(tempMPS,(),(1,2,3)),V
end

