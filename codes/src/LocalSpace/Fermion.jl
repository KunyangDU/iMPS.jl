module SpinlessFermion
using TensorKit

const d = 2
const PhySpace = (ℂ^d)'
const JWstring2 = TensorMap([-1 0; 0 1],PhySpace,PhySpace)
const FermionCreation2 = TensorMap([0 1;0 0],PhySpace,PhySpace)
const FermionAnnihilation2 = TensorMap([0 0;1 0],PhySpace,PhySpace)

const Z = let 
    JWstring2
end

const FFdag = let 
    FermionAnnihilation2,FermionCreation2
end

const FdagF = let 
    FermionCreation2,FermionAnnihilation2
end

const n = let 
    FermionCreation2*FermionAnnihilation2
end

const nn = let 
    n, n
end

end

module Spin2Fermion
using TensorKit

function diagm(dg::Vector{T}) where T
    L = length(dg)
    mat = zeros(T,L,L)
    for (dgi,dge) in enumerate(dg)
        mat[dgi,dgi] = dge
    end
    return mat
end

function diagm(pair::Pair{Int64, Vector{T}}) where T
    L = length(pair[2]) + abs(pair[1])
    mat = zeros(T,L,L)
    if pair[1] > 0
        for (ii,ie) in enumerate(pair[2])
            mat[ii,ii+pair[1]] = ie
        end
    elseif pair[1] < 0
        for (ii,ie) in enumerate(pair[2])
            mat[ii-pair[1],ii] = ie
        end
    else
        mat = diagm(pair[2])
    end
    
    return mat
end

const d = 4
const PhySpace = (ℂ^d)'

const JWstring4 = TensorMap(diagm([1,-1,-1,1]),PhySpace,PhySpace)
const FermionCreationUp4 = TensorMap(diagm(2 => [1,1]),PhySpace,PhySpace)
const FermionCreationDown4 = TensorMap(diagm(1 => [1,0,1]),PhySpace,PhySpace)

const Z = let 
    JWstring4
end

const FFdagUp = let 
    -FermionCreationUp4'*JWstring4,FermionCreationUp4
end

const FFdagDown = let 
    FermionCreationDown4',-JWstring4*FermionCreationDown4
end

const FdagFUp = let 
    FermionCreationUp4*JWstring4,FermionCreationUp4'
end

const FdagFDown = let 
    FermionCreationDown4,JWstring4*FermionCreationDown4'
end

const nup = let 
    FermionCreationUp4*FermionCreationUp4'
end

const ndown = let 
    FermionCreationDown4*FermionCreationDown4'
end

const n = let 
    nup + ndown
end

const nd = let 
    nup*ndown
end

const SS = let 
    (nup - ndown,nup - ndown)
end

end
