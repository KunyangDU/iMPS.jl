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