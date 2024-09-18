
function RandMPS(L::Int64;d::Int64=2)
    # suppose center at leftmost
    MPS = Vector{AbstractTensorMap}(undef,L)

    bond = ℂ^1
    phys = (ℂ^d)'

    MPS[1] = Tensor(rand(ComplexF64,1,d,1),bond' ⊗ phys' ⊗ bond') |> x -> permute(x,(),(1,2,3)) / norm(x)
    
    for i in 2:L
        MPS[i] = let
            TensorMap(rand(1,d,1),bond,phys ⊗ bond) |> x -> x / norm(x)
        end
    end
    
    return MPS
end

function ApproxReal(Qi::Number;tol::Float64=1e-5)
    imag(Qi) <= tol && return real(Qi)
    @error "not real"
end


function LocalMerge(ψ1::AbstractTensorMap{ComplexSpace,1,2},
    ψ2::AbstractTensorMap{ComplexSpace,0,3})

    @tensor ψm[-1,-2,-3,-4] ≔ ψ1[1,-1,-2] * ψ2[1,-3,-4]
    return permute(ψm,(),(1,2,3,4))
end

function LocalMerge(ψ1::AbstractTensorMap{ComplexSpace,0,3},
    ψ2::AbstractTensorMap{ComplexSpace,1,2})

    @tensor ψm[-1,-2,-3,-4] ≔ ψ1[-1,-2,1] * ψ2[1,-3,-4]
    return permute(ψm,(),(1,2,3,4))
end

function LocalMerge(
    EnvL::AbstractTensorMap{ComplexSpace,1,1},
    Opr1::AbstractTensorMap{ComplexSpace,1,3},
    Opr2::AbstractTensorMap{ComplexSpace,2,2},
    EnvR::AbstractTensorMap{ComplexSpace,1,1}
    )

    @tensor tempMPO[-1,-2,-3,-4,-5,-6] ≔ EnvL[1,-3]*Opr1[-2,1,-4,2]*Opr2[-1,2,-5,3]*EnvR[3,-6]
    return permute(tempMPO,(1,2),(3,4,5,6))
end

function LocalMerge(
    EnvL::AbstractTensorMap{ComplexSpace,1,1},
    Opr1::AbstractTensorMap{ComplexSpace,2,2},
    Opr2::AbstractTensorMap{ComplexSpace,1,3},
    EnvR::AbstractTensorMap{ComplexSpace,1,1}
    )

    @tensor tempMPO[-1,-2,-3,-4,-5,-6] ≔ EnvL[1,-3]*Opr1[2,-2,1,-4]*Opr2[-1,2,-5,3]*EnvR[3,-6]
    return permute(tempMPO,(1,2),(3,4,5,6))
end


function InnerProd(ψ₁::Vector,ψ₂::Vector)
    EnvL = InitialLeftEnv(;order=2)
    EnvR = RightEnv(ψ₁,ψ₂,1)
    innerprod = @tensor EnvL[1,3]*ψ₁[1][1,5,2]*ψ₂[1]'[3,5,4]*EnvR[2,4]
    return innerprod[1]
end

function vrange(beginvec::Union{Vector,Tuple},endvec::Union{Vector,Tuple};step::Int64 = 100)
    return hcat([collect(beginvec .+ (endvec .- beginvec) .* t)  for t in range(0,1,step)]...)
end

function vrange(ipath::Matrix;eachstep::Number = 100)
    finalpath = ipath[:,1]

    for ii in 1:size(ipath)[2]-1
        finalpath = hcat(finalpath,vrange(ipath[:,ii],ipath[:,ii+1];step=eachstep+1)[:,2:end])
    end

    return finalpath
end

function pathlength(finalpath::Matrix)
    return cumsum(norm.(eachcol(hcat([0.0,0.0],diff(finalpath,dims = 2)))))
end

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

function kdivide(kvecpath::Matrix,groupn::Int64)
    nperg = div(size(kvecpath)[2] - 1,groupn)
    kg = []
    for ii in 1:groupn-1
        push!(kg,kvecpath[:,(ii-1)*nperg .+ (1:nperg)])
    end
    push!(kg,kvecpath[:,end-nperg:end])

    return kg
end

function kdivide(kr::Vector,groupn::Int64)
    nperg = div(length(kr) - 1,groupn)
    kg = []
    for ii in 1:groupn-1
        push!(kg,kr[(ii-1)*nperg .+ (1:nperg)])
    end
    push!(kg,kr[end-nperg:end])
    return kg
end


function canonicalize(Opr::Vector{AbstractTensorMap{ComplexSpace,2,2}})

    Opr = convert(Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}},Opr)
    Opr[end] = Move(Opr[end],InitialEnv(Vector{AbstractTensorMap{ComplexSpace,1,3}}(undef,1)))[1]
    for i in length(Opr):-1:2
        Opr[i-1:i] = Move(Opr[i-1:i]...)
    end

    return Opr
end

