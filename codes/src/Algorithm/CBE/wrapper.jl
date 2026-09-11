function CBE!(env::Environment, alg::CBEalgo{dynamicSVD}, info::CBEinfo{Dir};kwargs...) where Dir
    to = TimerOutput()

    Dl,Dr = _cbe_maxdim(env,alg,info)
    Dc = _cbe_currentdim(env,alg,info)

    if Dl > alg.D && Dr > alg.D
        # full but needs trunc
        @timeit to "rand SVD" localto = CBE!(env,CBEalgo(alg,randSVD()),info)
        merge!(to,localto,tree_point = ["rand SVD"])
        return to
    elseif Dl > Dc && Dr > Dc
        # not full yet
        @timeit to "full SVD" localto = CBE!(env,CBEalgo(alg,fullSVD()),info)
        merge!(to,localto,tree_point = ["full SVD"])
        return to
    else
        # full but no trunc is need
        return to
    end
end

CBE!(env::CBEenvironment, alg::CBEalgo{randSVD}, info::CBEinfo;kwargs...) = randSVD!(env,alg,info)
CBE!(env::CBEenvironment, alg::CBEalgo{fullSVD}, info::CBEinfo;kwargs...) = fullSVD!(env,alg,info)

