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



