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
