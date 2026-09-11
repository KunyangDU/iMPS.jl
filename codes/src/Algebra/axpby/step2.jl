
function axpby!(α::Number, Envx::Environment{2}, β::Number, Envy::Environment{2}, Alg::Algebraalgo{DoubleSite}, sweepinfo::Algebrasweepinfo{L2R};kwargs...)
    localto = TimerOutput()
    L = length(Envx.layer[1])
    for site in 1:L-1
        Alg.verbose && (time₀ = time())
        localinfo = Algebrasiteinfo()
        x₀ = deepcopy(composite(Envx.layer[1][site:site+1]...))
        @assert (x2 = norm(x₀)^2) ≠ 0
        @timeit localto "action_x" tx = actionb(proj2(Envx.envs[site], nothing, nothing, Envx.envs[site+2]), Envx.layer[2][site:site+1]...)
        @timeit localto "action_y" ty = actionb(proj2(Envy.envs[site], nothing, nothing, Envy.envs[site+2]), Envy.layer[2][site:site+1]...)
        # 用 add!! 而非就地 axpby!：tx（来自 Envx，复）与 ty（来自 Envy，实）标量类型可不同（如复 H 的有限温），
        # 就地写 ty.A 会 InexactError（与 _sparse_actionb_sum 的提升同理）。
        @timeit localto "SVD" tl, tr, localinfo.err, localinfo.bond = tsvd(add!!(tx, ty, β, α)'; direction=:right,trunc = Alg.trunc)
        @timeit localto "push right" map([Envx,Envy]) do Env
            Env.layer[1][site:site+1] = tl, tr
            canonicalize!!(Env.layer[1], site+1)
            canonicalize!(Env, site+1)
        end
        x = composite(Envx.layer[1][site:site+1]...)
        localinfo.err = norm(x-x₀)^2/x2

        merge!(sweepinfo,localinfo)
        Alg.verbose && vbshow(site, time₀, localinfo, Alg)
    end
    return localto
end

function axpby!(α::Number, Envx::Environment{2}, β::Number, Envy::Environment{2}, Alg::Algebraalgo{DoubleSite}, sweepinfo::Algebrasweepinfo{R2L};kwargs...)
    localto = TimerOutput()
    L = length(Envx.layer[1])
    for site in L:-1:2
        Alg.verbose && (time₀ = time())
        localinfo = Algebrasiteinfo()
        x₀ = deepcopy(composite(Envx.layer[1][site-1:site]...))
        @assert (x2 = norm(x₀)^2) ≠ 0
        @timeit localto "action" tx = actionb(proj2(Envx.envs[site-1], nothing, nothing, Envx.envs[site+1]), Envx.layer[2][site-1:site]...)
        @timeit localto "action" ty = actionb(proj2(Envy.envs[site-1], nothing, nothing, Envy.envs[site+1]), Envy.layer[2][site-1:site]...)
        # 同 L2R：add!! 做类型提升，避免实/复混用时的 InexactError。
        @timeit localto "SVD" tl, tr, localinfo.err, localinfo.bond = tsvd(add!!(tx, ty, β, α)'; direction=:left,trunc = Alg.trunc)
        @timeit localto "push left" map([Envx,Envy]) do Env
            Env.layer[1][site-1:site] = tl, tr
            canonicalize!!(Env.layer[1], site-1)
            canonicalize!(Env, site-1)
        end
        x = composite(Envx.layer[1][site-1:site]...)
        localinfo.err = norm(x-x₀)^2/x2

        merge!(sweepinfo,localinfo)
        Alg.verbose && vbshow(site, time₀, localinfo, Alg)
    end
    return localto
end