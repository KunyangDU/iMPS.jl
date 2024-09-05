

function update1(ψ::Vector,H::Vector,site::Int64,level::Int64,direction::String)

    reduH = reduHam(ψ,H,site)
    Eg,Ev = groundEig(reduH,level)
    MPSs = orientSVD(ψ,Ev',site,direction)

#=     if site == 1
        #= HR = RightEnv(ψ,H,site)
        @tensor reduH[-1,-2,-3,-4] ≔ H[site][-1,1,-3]*HR[-2,1,-4]
        reduH = permute(reduH,(1,2),(3,4)) =#

        Eg,Ev = groundEig(reduH,level)

        #= U,S,V = tsvd(Ev',(2,),(1,))
        @tensor tempMPS[-1,-2,-3] ≔ permute(U*S,(),(1,2))[1,-1]*ψ[2][1,-2,-3]
        nextMPS = permute(tempMPS,(),(1,2,3))
        thisMPS = V =#
    elseif site == length(ψ)
        #= HL = LeftEnv(ψ,H,site)
        showdomain(H[L])
        showdomain(HL)
        showdomain(H[site])
        @tensor reduH[-1,-2,-3,-4] ≔ HL[]*H[site][-1,1,-3] =#

        Eg,Ev = groundEig(reduH,level)
    else
        #= HR = RightEnv(ψ,H,site)
        HL = LeftEnv(ψ,H,site)
        @tensor reduH[-1,-2,-3,-4,-5,-6] ≔ HL[-1,1,-4]*H[site][-2,2,-5,1]*HR[-3,2,-6]
        reduH = permute(reduH,(1,2,3),(4,5,6)) =#
        
        Eg,Ev = groundEig(reduH,level)

        #= U,S,V = tsvd(Ev',(3,),(1,2))
        thisMPS = V
        if site == L-1
            @tensor tempMPS[-1,-2] ≔ permute(U*S,(),(1,2))[-1,1]*ψ[site+1][1,-2]
            nextMPS = permute(tempMPS,(),(1,2))
        else
            @tensor tempMPS[-1,-2,-3] ≔ permute(U*S,(),(1,2))[-1,1]*ψ[site+1][1,-2,-3]
            nextMPS = permute(tempMPS,(),(1,2,3))
        end =#
    end =#

    return Eg, MPSs
end



function sweep1(ψ::Vector,H::Vector,Nsweep::Int64,LanczosLevel::Int64,D_MPS::Int64)
    lsE = []

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
        push!(lsE,Eg)
    end

    return ψ,lsE

end




