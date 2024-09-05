

function showdomain(M::AbstractTensorMap)
    @show codomain(M),domain(M)
end

function showQuantSweep(lsQ::Vector;name::String = "Quantity")
    for (iq,q) in enumerate(lsQ)
        println("$name\t$iq\t$q")
    end
end

