
function RandMPS(L::Int64;d::Int64=2)
    # suppose center at leftmost
    MPS = Vector{Union{AbstractTensorMap{ComplexSpace,0,3},AbstractTensorMap{ComplexSpace,1,2}}}(undef,L)

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
    EnvL::AbstractTensorMap{ComplexSpace,2,1},
    Opr11::AbstractTensorMap{ComplexSpace,1,3},
    Opr12::AbstractTensorMap{ComplexSpace,2,2},
    Opr21::AbstractTensorMap{ComplexSpace,1,3},
    Opr22::AbstractTensorMap{ComplexSpace,2,2},
    EnvR::AbstractTensorMap{ComplexSpace,2,1}
    )

    @tensor tempMPO[-1,-2,-3,-4,-5,-6] ≔ EnvL[1,2,-3]*Opr11[-2,1,3,4]*Opr21[3,2,-4,5]*Opr12[-1,4,6,7]*Opr22[6,5,-5,8]*EnvR[7,8,-6]
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

function LocalMerge(
    EnvL::AbstractTensorMap{ComplexSpace,2,1},
    Opr11::AbstractTensorMap{ComplexSpace,2,2},
    Opr12::AbstractTensorMap{ComplexSpace,1,3},
    Opr21::AbstractTensorMap{ComplexSpace,2,2},
    Opr22::AbstractTensorMap{ComplexSpace,1,3},
    EnvR::AbstractTensorMap{ComplexSpace,2,1}
    )

    @tensor tempMPO[-1,-2,-3,-4,-5,-6] ≔ EnvL[7,8,-3]*Opr11[4,-2,7,6]*Opr21[5,6,8,-4]*Opr12[-1,4,3,1]*Opr22[3,5,-5,2]*EnvR[1,2,-6]
    return permute(tempMPO,(1,2),(3,4,5,6))
end

function LocalMerge(
    EnvL::AbstractTensorMap{ComplexSpace,1,1},
    ψ1::AbstractTensorMap{ComplexSpace,0,3},
    ψ2::AbstractTensorMap{ComplexSpace,1,2},
    EnvR::AbstractTensorMap{ComplexSpace,1,1}
    )

    @tensor tempMPS[-1,-2,-3,-4] ≔ EnvL[1,-1]*ψ1[1,-2,2]*ψ2[2,-3,3]*EnvR[3,-4]
    return permute(tempMPS,(),(1,2,3,4))
end

function LocalMerge(
    EnvL::AbstractTensorMap{ComplexSpace,1,1},
    ψ1::AbstractTensorMap{ComplexSpace,1,2},
    ψ2::AbstractTensorMap{ComplexSpace,0,3},
    EnvR::AbstractTensorMap{ComplexSpace,1,1}
    )

    @tensor tempMPS[-1,-2,-3,-4] ≔ EnvL[3,-1]*ψ1[2,3,-2]*ψ2[2,-3,1]*EnvR[1,-4]
    return permute(tempMPS,(),(1,2,3,4))
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
    #= 
    canonicalize!
    从SETTN开始所有的MPO都是正则化的，包括RandMPO
    其它零温计算的MPO都不正则化
    =#
    Opr = convert(Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}},Opr)
    Opr[end] = Move(Opr[end],InitialEnv(Vector{AbstractTensorMap{ComplexSpace,1,3}}(undef,1)))[1]
    for i in length(Opr):-1:2
        Opr[i-1:i] = Move(Opr[i-1:i]...)
    end

    return Opr
end


function IdentityMPO(L::Int64,d::Int64)
    
    MPO = Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}}(undef,L)

    idt = ℂ^1
    phys = (ℂ^d)'

    MPO[1] = TensorMap(diagm(ones(d)),phys,idt ⊗ phys ⊗ idt)
    
    for i in 2:L
        MPO[i] = let
            TensorMap(diagm(ones(d)),phys ⊗ idt,phys ⊗ idt)
        end
    end
    println("MPO constructed")

    return MPO
end


function Trace(Opr1::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}})
    
    EnvR = RightEnv(Opr1,Opr1,1)
    EnvL = LeftEnv(Opr1,Opr1,1)

    tr = @tensor EnvL[1,2]*Opr1[1][4,1,3,5]*Opr1[1]'[2,3,6,4]*EnvR[5,6]

    return ApproxReal(tr[1])
end

function Trace(
    Opr1::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}},
    Opr2::Vector{Union{AbstractTensorMap{ComplexSpace,2,2},AbstractTensorMap{ComplexSpace,1,3}}};
    Z::Number=1
    )

    EnvR = RightEnv(Opr1,Opr2,Opr1,1)
    EnvL = LeftEnv(Opr1,Opr2,Opr1,1)

    tr = @tensor EnvL[1,2,4]*Opr1[1][6,1,3,7]*Opr2[1][3,2,5,8]*Opr1[1]'[4,5,9,6]*EnvR[7,8,9]

    return ApproxReal(tr[1]) / Z
end


function WeightSum(MPSs::Vector{AbstractTensorMap},weights::Vector,D_MPS::Int64)
    d = dims(domain(MPSs[1]))[2]
    L = length(MPSs[1])

    MPS = VariPlusMPS(vcat([weights[1]],ones(L-1)) .* MPSs[1], vcat([weights[2]],ones(L-1)) .* MPSs[2],d,D_MPS)

    for i in eachindex(MPSs)[3:end]
        MPS = VariPlusMPS(MPS,vcat([weights[i]],ones(L-1)) .* MPSs[i],d,D_MPS)
    end
    
    return MPS
end
