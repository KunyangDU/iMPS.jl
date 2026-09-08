

function axpby!(α::Number, x::DenseMPO{L}, β::Number, y::DenseMPO{L}, Alg::Algebraalgo;kwargs...) where L
    
    y′ = y'
    to = TimerOutput()
    __init_io__()
    @timeit to "initialize XY Env" begin
        Envx = Environment([y,RefMPO(x,adjoint)];isdisk=Alg.isdisk)
        Envy = Environment([y,y′];isdisk=Alg.isdisk)
        initialize!(Envx)
        initialize!(Envy)
    end

    info = Algebrainfo()
    try
        while info.n ≤ Alg.N
            localto = TimerOutput()

            l2rinfo = Algebrasweepinfo(L2R())
            mto = axpby!(α,Envx,β,Envy,Alg,l2rinfo)

            show(mto;title = ">>> axpby! >>>")
            print("\n")
            show(l2rinfo)
            flush(stdout)

            merge!(localto,mto)
            merge!(info,l2rinfo)

            r2linfo = Algebrasweepinfo(R2L())
            mto = axpby!(α,Envx,β,Envy,Alg,r2linfo)

            show(mto;title = "<<< axpby! <<<")
            print("\n")
            show(r2linfo)
            flush(stdout)

            merge!(localto,mto)
            merge!(info,r2linfo)

            _merge_io!(localto)
            merge!(to,localto)

            info.err < Alg.tol && break
        end
        return y, to, info
    finally
        Alg.isdisk && (cleanup!(Envx); cleanup!(Envy); cleanup!(y′))
    end
end

axpy!(α::Number, x::DenseMPO, y::DenseMPO;kwargs...) = axpby!(α,x,1,y;kwargs...)
axpy!(α::Number, x::DenseMPO, y::DenseMPO, algo::Algebraalgo;kwargs...) = axpby!(α,x,1,y,algo;kwargs...)

