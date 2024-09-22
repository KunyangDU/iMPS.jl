function RightMove(nextψ::AbstractTensorMap{ComplexSpace, 1, 2},thisψ::AbstractTensorMap{ComplexSpace, 0, 3},D_MPS::Int64)
    MPSs,truncerr = RightSVD(thisψ,D_MPS)
    V,Σ = MPSs
    return RightMerge(Σ,V,nextψ),truncerr
end

function Move(
    Opri1::Union{AbstractTensorMap{ComplexSpace, 1, 3},AbstractTensorMap{ComplexSpace, 2, 2}},
    Opri2::Union{AbstractTensorMap{ComplexSpace, 1, 3},AbstractTensorMap{ComplexSpace, 2, 2}}
    )
    Opris = mySVD(Opri1,Opri2)
    return Merge(Opri1,Opri2,Opris...)
end

function Move(
    ψ1::Union{AbstractTensorMap{ComplexSpace, 0, 3},AbstractTensorMap{ComplexSpace, 1, 2}},
    ψ2::Union{AbstractTensorMap{ComplexSpace, 0, 3},AbstractTensorMap{ComplexSpace, 1, 2}},
    )
    MPSs = mySVD(ψ1,ψ2)
    return Merge(ψ1,ψ2,MPSs...)
end

function LeftMove(nextψ::AbstractTensorMap{ComplexSpace, 1, 2},thisψ::AbstractTensorMap{ComplexSpace, 0, 3},D_MPS::Int64)
    MPSs,truncerr = LeftSVD(thisψ,D_MPS)
    Σ,V = MPSs
    return LeftMerge(Σ,V,nextψ),truncerr
end

# SVD

function RightSVD(ψi::AbstractTensorMap{ComplexSpace, 0, 3},D_MPS::Int64)
    U,S,V,truncerr = tsvd(ψi,(3,),(1,2);trunc = truncdim(D_MPS))
    return [V,permute(U*S,(),(1,2))],truncerr
end

function RightSVD(ψi::AbstractTensorMap{ComplexSpace, 0, 4},D_MPS::Int64)
    U,S,V,truncerr = tsvd(ψi,(3,4),(1,2);trunc = truncdim(D_MPS))
    return [V,permute(U*S,(),(3,1,2))],truncerr
end

function LeftSVD(ψi::AbstractTensorMap{ComplexSpace,0,3},D_MPS::Int64)
    U,S,V,truncerr = tsvd(ψi,(1,),(2,3);trunc = truncdim(D_MPS))
    return [permute(U*S,(),(1,2)),permute(V,(1,),(2,3))],truncerr
end

function LeftSVD(ψi::AbstractTensorMap{ComplexSpace,0,4},D_MPS::Int64)
    U,S,V,truncerr = tsvd(ψi,(1,2),(3,4);trunc = truncdim(D_MPS))
    return [permute(U*S,(),(1,2,3)),permute(V,(1,),(2,3))],truncerr
end

function mySVD(
    Opri1::AbstractTensorMap{ComplexSpace, 1, 3},
    Opri2::AbstractTensorMap{ComplexSpace, 2, 2}
    )
    U,S,V = tsvd(Opri1,(4,),(1,2,3))
    return [permute(V,(1,2),(3,4)),permute(U*S,(),(1,2))]
end

function mySVD(
    Opri1::AbstractTensorMap{ComplexSpace, 2, 2},
    Opri2::AbstractTensorMap{ComplexSpace, 1, 3}
    )
    U,S,V = tsvd(Opri2,(2,),(3,4,1))
    return [permute(U*S,(),(1,2)),permute(V,(4,1),(2,3))]
end

# 如果2个非正则MPO，则默认向左SVD
function mySVD(
    Opri1::AbstractTensorMap{ComplexSpace, 2, 2},
    Opri2::AbstractTensorMap{ComplexSpace, 2, 2}
    )
    U,S,V = tsvd(Opri2,(3,),(4,1,2))
    showdomain(V)
    return [permute(U*S,(),(1,2)),permute(V,(4,1),(2,3))]
end

function mySVD(Opri::AbstractTensorMap{ComplexSpace, 2, 4},
    direction::String, D_MPO::Int64)
    if direction == "right"
        U,S,V,truncerr = tsvd(Opri,(5,6,1),(2,3,4);trunc = truncdim(D_MPO))
        return [permute(V,(1,2),(3,4)),permute(U*S,(3,),(4,1,2))],truncerr
    elseif direction == "left"
        U,S,V,truncerr = tsvd(Opri,(2,3,4),(5,6,1);trunc = truncdim(D_MPO))
        return [permute(U*S,(1,),(2,3,4)),permute(V,(4,1),(2,3))],truncerr
    else
        @error "direction does not exist!"
    end
end


function mySVD(
    ψ1::AbstractTensorMap{ComplexSpace, 0, 3},
    ψ2::AbstractTensorMap{ComplexSpace, 1, 2}
    )
    U,S,V = tsvd(ψ1,(3,),(1,2))
    return [V,permute(U*S,(),(1,2))]
end

function mySVD(
    ψ1::AbstractTensorMap{ComplexSpace, 1, 2},
    ψ2::AbstractTensorMap{ComplexSpace, 0, 3}
    )
    U,S,V = tsvd(ψ2,(1,),(2,3))
    return [permute(U*S,(),(1,2)),V]
end

function mySVD(ψi::AbstractTensorMap{ComplexSpace, 0, 4},
    direction::String, D_MPS::Int64)
    if direction == "right"
        U,S,V,truncerr = tsvd(ψi,(3,4),(1,2);trunc = truncdim(D_MPS))
        return [V,permute(U*S,(),(3,1,2))],truncerr
    elseif direction == "left"
        U,S,V,truncerr = tsvd(ψi,(1,2),(3,4);trunc = truncdim(D_MPS))
        return [permute(U*S,(),(1,2,3)),V],truncerr
    else
        @error "direction does not exist!"
    end
end

# MERGE

function RightMerge(Σ::AbstractTensorMap{ComplexSpace, 0, 2},V::AbstractTensorMap{ComplexSpace, 1, 2},nextψ::AbstractTensorMap{ComplexSpace, 1, 2})
    @tensor tempMPS[-1,-2,-3] ≔ Σ[1,-1]*nextψ[1,-2,-3]
    return [V,permute(tempMPS,(),(1,2,3))]
end

function LeftMerge(Σ::AbstractTensorMap{ComplexSpace, 0, 2},V::AbstractTensorMap{ComplexSpace, 1, 2},nextψ::AbstractTensorMap{ComplexSpace, 1, 2})
    @tensor tempMPS[-1,-2,-3] ≔ Σ[1,-3]*nextψ[1,-1,-2]
    return [permute(tempMPS,(),(1,2,3)),V]
end

function Merge(
    Opri1::AbstractTensorMap{ComplexSpace, 1, 3},
    Opri2::AbstractTensorMap{ComplexSpace, 2, 2},
    svdOpri1::AbstractTensorMap{ComplexSpace, 2, 2},
    svdOpri2::AbstractTensorMap{ComplexSpace, 0, 2},
    )
    @tensor tempMPO[-1,-2,-3,-4] ≔ svdOpri2[1,-2]*Opri2[-1,1,-3,-4]

    return [svdOpri1,permute(tempMPO,(1,),(2,3,4))]
end

function Merge(
    Opri1::AbstractTensorMap{ComplexSpace, 2, 2},
    Opri2::AbstractTensorMap{ComplexSpace, 1, 3},
    svdOpri1::AbstractTensorMap{ComplexSpace, 0, 2},
    svdOpri2::AbstractTensorMap{ComplexSpace, 2, 2},
    )
    @tensor tempMPO[-1,-2,-3,-4] ≔ Opri1[1,-1,-2,-3]*svdOpri1[1,-4]

    return [permute(tempMPO,(1,),(2,3,4)),svdOpri2]
end

function Merge(
    Opri1::AbstractTensorMap{ComplexSpace, 2, 2},
    Opri2::AbstractTensorMap{ComplexSpace, 2, 2},
    svdOpri1::AbstractTensorMap{ComplexSpace, 0, 2},
    svdOpri2::AbstractTensorMap{ComplexSpace, 2, 2},
    )
    @tensor tempMPO[-1,-2,-3,-4] ≔ Opri1[1,-1,-2,-3]*svdOpri1[1,-4]

    return [permute(tempMPO,(1,),(2,3,4)),svdOpri2]
end

function Merge(
    ψ1::AbstractTensorMap{ComplexSpace, 0, 3},
    ψ2::AbstractTensorMap{ComplexSpace, 1, 2},
    svdψ1::AbstractTensorMap{ComplexSpace, 1, 2},
    svdψ2::AbstractTensorMap{ComplexSpace, 0, 2},
    )
    @tensor tempMPS[-1,-2,-3] ≔ svdψ2[1,-1]*ψ2[1,-2,-3]

    return [svdψ1,permute(tempMPS,(),(1,2,3))]
end

function Merge(
    ψ1::AbstractTensorMap{ComplexSpace, 1, 2},
    ψ2::AbstractTensorMap{ComplexSpace, 0, 3},
    svdψ1::AbstractTensorMap{ComplexSpace, 0, 2},
    svdψ2::AbstractTensorMap{ComplexSpace, 1, 2},
    )
    @tensor tempMPS[-1,-2,-3] ≔ svdψ1[1,-3]*ψ1[1,-1,-2]

    return [permute(tempMPS,(),(1,2,3)),svdψ2]
end
