function fullSVD!(env::CBEenvironment, alg::CBEalgo,info::CBEinfo{L2R})
    localto = TimerOutput()
    @timeit localto "composite" x = composite(env.tL₀,env.tR₀)
    @timeit localto "rightorth" env.tL,env.tR,~,info.bond = tsvd(x;direction = :left,trunc = truncdim(alg.D))
    return localto
end

function fullSVD!(env::CBEenvironment,alg::CBEalgo,info::CBEinfo{R2L})
    localto = TimerOutput()
    @timeit localto "composite" x = composite(env.tL₀,env.tR₀)
    @timeit localto "leftorth" env.tL,env.tR,~,info.bond  = tsvd(x;direction = :right,trunc = truncdim(alg.D))
    return localto
end

function randSVD!(env::CBEenvironment, alg::CBEalgo,info::CBEinfo{L2R})
    localto = TimerOutput()

    Ω = _cbetensor(randn,env.tR₀,alg.D,L2R())

    @timeit localto "splice Ω" R_trunc = splice(env.Rorth,Ω)
    @timeit localto "contract_LO*Rt" Q = contract(env.Lorth, R_trunc, env.lm)
    @timeit localto "leftorth" Q,~ = leftorth(Q)
    @timeit localto "splice Q'" L_trunc = splice(env.Lorth,Q')
    @timeit localto "contract_Lt*RO" obj = contract(L_trunc, env.Rorth, env.lm)
    @timeit localto "SVD" ~,tR′,info.err,info.bond = tsvd(obj';direction = :left,trunc = truncdim(alg.D - dims(env.tL₀)[2][1]) & truncbelow(alg.tol))

    @timeit localto "oplus" begin 
        tL′, tR′ = _rexpand(env.tL₀, tR′)
        env.tL = _roplus(env.tL₀, tL′)
        env.tR = _loplus(env.tR₀, tR′)
    end

    # @timeit localto "check" begin
    #     x₀ = composite(env.tL₀,env.tR₀)
    #     x = composite(env.tL,env.tR)
    #     @show norm(x - x₀)^2
    #     # @show norm(x),norm(x₀)
    # end

    return localto
end

function randSVD!(env::CBEenvironment,alg::CBEalgo,info::CBEinfo{R2L})
    localto = TimerOutput()

    Ω = _cbetensor(randn,env.tL₀,alg.D,R2L())

    @timeit localto "splice Ω" L_trunc = splice(env.Lorth,Ω)
    @timeit localto "contract_Lt*RO" Q = contract(L_trunc, env.Rorth, env.lm)
    @timeit localto "rightorth" ~,Q = rightorth(Q)
    @timeit localto "splice Q'" R_trunc = splice(env.Rorth,Q')
    @timeit localto "contract_LO*Rt" obj = contract(env.Lorth, R_trunc, env.lm)
    @timeit localto "SVD" tL′,~,info.err,info.bond = tsvd(obj';direction = :right,trunc = truncdim(alg.D - dims(env.tL₀)[2][1]) & truncbelow(alg.tol))

    @timeit localto "oplus" begin 
        tL′, tR′ = _lexpand(tL′, env.tR₀)
        env.tL = _roplus(env.tL₀, tL′)
        env.tR = _loplus(env.tR₀, tR′)
    end

    # @timeit localto "check" begin
    #     x₀ = composite(env.tL₀,env.tR₀)
    #     x = composite(env.tL,env.tR)
    #     @show norm(x - x₀)^2
    #     # @show norm(x),norm(x₀)
    # end

    return localto
end


