function showdomain(M::AbstractTensorMap)
    @show codomain(M),domain(M)
end

function showQuantSweep(lsQ::Vector;name::String = "Quantity")
    for (iq,q) in enumerate(lsQ)
        println("$name\t$iq\t$q")
    end
end

function showBlockMPO(Block::AbstractTensorMap{ComplexSpace,2,2})
    # where block of MPO is shown in block matrix, and the outer index is bond
    # i.e.
    #  I  I₀
    #  σₓ I 
    ChosenBlock = permute(Block,(2,4),(3,1))
    @show ChosenBlock
end

function showMatrixElem(M::Matrix)
    row,col = size(M)
    
    fig = Figure()
    ax = Axis(fig[1,1],
    xticks = 1:row,
    yticks = 1:col)
    
    heatmap!(ax,M'[:,end:-1:1])
    
    display(fig)    
end

function checkOrth(M::AbstractTensorMap{ComplexSpace,1,2})
    @tensor test[-1,-2] ≔ M[-1,1,2]*M'[1,2,-2]
    return test
end

function irange(start::Int64,stop::Int64)
    if start > stop
        return start:-1:stop 
    else
        return start:stop
    end
end

function FindMaxDist(pairs::Vector)
    return maximum(abs(i - j) for (i, j) in pairs)
end

function FindPair(pairs::Vector,site::Number,target::Number)
    return filter(x -> x[site] == target, pairs)
end
