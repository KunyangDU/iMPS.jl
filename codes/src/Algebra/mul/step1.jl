
function mul!(EnvAB::Environment{3}, α::Number, Alg::Algebraalgo{SingleSite,alg}, sweepinfo::Algebrasweepinfo{L2R}; kwargs...) where alg
    localto = TimerOutput()
    L = length(EnvAB.layer[1])
    for site in 1:L-1
        Alg.verbose && (time₀ = time())
        localinfo = Algebrasiteinfo()
        x₀ = composite((EnvAB.layer[1][site:site+1])...)
        @assert (x2 = norm(x₀)^2) ≠ 0

        if alg <: CBEalgo 
            cbeinfo = CBEinfo(L2R())
            @timeit localto "CBE_AB" cbetoAB = CBE!(EnvAB, Alg.alg, cbeinfo)
            # merge!(localinfo,cbeinfo)
            merge!(localto,cbetoAB,tree_point = ["CBE_AB"])
        end

        @timeit localto "projection" projH = proj1(EnvAB,site)
        @timeit localto "action" t = actionb(projH,EnvAB.layer[3][site])
        Norm = normalize!(rmul!(t, α))
        @timeit localto "svd" EnvAB.layer[1][site],tr,localinfo.truncerr,localinfo.bond = tsvd(t'; direction=:right,trunc = Alg.trunc)
        normalize!(tr)
        EnvAB.layer[1][site+1] = splice(rmul!(tr, Norm),EnvAB.layer[1][site+1])
        canonicalize!!(EnvAB.layer[1],site+1)
        @timeit localto "push right" canonicalize!(EnvAB, site+1)

        x = composite((EnvAB.layer[1][site:site+1])...)
        localinfo.err = norm(x-x₀)^2/x2

        merge!(sweepinfo,localinfo)
        Alg.verbose && vbshow(site, time₀, localinfo, Alg)
    end

    return localto
end

function mul!(EnvAB::Environment{3}, α::Number, Alg::Algebraalgo{SingleSite,alg}, sweepinfo::Algebrasweepinfo{R2L}; kwargs...) where alg
    localto = TimerOutput()
    L = length(EnvAB.layer[1])
    for site in L:-1:2
        Alg.verbose && (time₀ = time())
        localinfo = Algebrasiteinfo()
        x₀ = composite((EnvAB.layer[1][site-1:site])...)
        @assert (x2 = norm(x₀)^2) ≠ 0

        if alg <: CBEalgo 
            cbeinfo = CBEinfo(R2L())
            @timeit localto "CBE_AB" cbetoAB = CBE!(EnvAB, Alg.alg, cbeinfo)
            # merge!(localinfo,cbeinfo)
            merge!(localto,cbetoAB,tree_point = ["CBE_AB"])
        end

        @timeit localto "projection" projH = proj1(EnvAB,site)
        @timeit localto "action" t = actionb(projH,EnvAB.layer[3][site])
        Norm = normalize!(rmul!(t, α))
        @timeit localto "svd" tl,EnvAB.layer[1][site],localinfo.truncerr,localinfo.bond = tsvd(t'; direction=:left,trunc = Alg.trunc)
        normalize!(tl)
        EnvAB.layer[1][site-1] = splice(EnvAB.layer[1][site-1],rmul!(tl, Norm))
        canonicalize!!(EnvAB.layer[1],site-1)
        @timeit localto "push left" canonicalize!(EnvAB,site-1)

        x = composite((EnvAB.layer[1][site-1:site])...)
        localinfo.err = norm(x-x₀)^2/x2
        merge!(sweepinfo,localinfo)
        Alg.verbose && vbshow(site, time₀, localinfo, Alg)
    end

    return localto
end




