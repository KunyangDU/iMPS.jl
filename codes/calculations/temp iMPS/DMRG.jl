
function sweepDMRG1(ψ::Vector,H::Vector,Nsweep::Int64,LanczosLevel::Int64,D_MPS::Int64)
    lsE = Vector{Float64}(undef,Nsweep)

    for i in 1:Nsweep
        println("sweep $i begin")
        start_time = time()
        Eg = 0
        env = Vector{AbstractTensorMap}(underf,2)
        for iL in 1:L-1
            Eg, ψ[iL:iL+1] = updateDMRG1(ψ,H,iL,LanczosLevel,"right",D_MPS)
            env = Vector{AbstractTensorMap}(underf,2)
        end
        println("right sweep finished")
        for iL in L:-1:2
            Eg, ψ[iL-1:iL] = updateDMRG1(ψ,H,iL,LanczosLevel,"left",D_MPS)
        end
        println("left sweep finished")
        println("sweep $i finished, Eg = $Eg, time consumed $(round(time()-start_time;digits=2))")
        lsE[i] = Eg
    end

    return ψ,lsE

end