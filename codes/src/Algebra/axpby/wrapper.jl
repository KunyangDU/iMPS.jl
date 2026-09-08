
"""
axpy!(α, x, y) -> y
Overwrite y with x * α + y and return y. If x and y have the same axes, it's equivalent with y .+= x .* a.
# kwargs
D_MPO: MPO bond dimension. Default is the maximum D of x, y.
Nsweep: times of variational calculation (sweep). Default is 2.

"""
function axpby!(α::Number, x::DenseMPO{L}, β::Number, y::DenseMPO{L};kwargs...) where L
    trunc = get(kwargs,:trunc,notrunc())
    N  = get(kwargs,:N,3)
    tol = get(kwargs,:tol,1e-8)
    verbose = get(kwargs,:verbose,false)
    isdisk = get(kwargs,:isdisk,IS_DISK[])
    algo = Algebraalgo(DoubleSite(),NoAlgorithm(),trunc,N,tol,verbose,isdisk)
    return axpby!(α,x,β,y,algo;kwargs...)
end