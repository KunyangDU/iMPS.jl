
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