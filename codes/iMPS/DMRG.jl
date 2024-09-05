

function update1(ψ::Vector,H::Vector,site::Int64,level::Int64,direction::String)

    reduH = reduHam(ψ,H,site)
    Eg,Ev = groundEig(reduH,level)
    MPSs = orientSVD(ψ,Ev',site,direction)

    return Eg, MPSs
end


function sweep1(ψ::Vector,H::Vector,Nsweep::Int64,LanczosLevel::Int64,D_MPS::Int64)
    lsE = Vector{Float64}(undef,Nsweep)

    for i in 1:Nsweep
        println("sweep $i begin")
        start_time = time()
        Eg = 0
        for iL in 1:L-1
            Eg, ψ[iL:iL+1] = update1(ψ,H,iL,LanczosLevel,"right")
        end
        println("right sweep finished")
        for iL in L:-1:2
            Eg, ψ[iL-1:iL] = update1(ψ,H,iL,LanczosLevel,"left")
        end
        println("left sweep finished")
        println("sweep $i finished, Eg = $Eg, time consumed $(round(time()-start_time;digits=2))")
        lsE[i] = Eg
    end

    return ψ,lsE

end




