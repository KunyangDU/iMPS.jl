function ChainBand(k::Number;t::Number = 1,a::Number=1)
    return 2*t*cos(k*a)
end

function SquaBand(k::Vector;t::Number = 1,a::Number=1)
    E = 0
    for ki in k 
        E += 2*t*cos(ki*a)
    end
    return E
end

function SquaBand(k::Number;t::Number = 1,a::Number=1)
    return 2*t*cos(k*a)
end


function TheoOccupCurve(;t::Number = 1,a::Number = 1,n::Int64 = 500)

    lskx = range(-pi/a,pi/a,n)
    lsky = range(-pi/a,pi/a,n)

    dkx = lskx[2]-lskx[1]
    dky = lsky[2]-lsky[1]

    E = zeros(n,n)

    for (i,kx) in enumerate(lskx),(j,ky) in enumerate(lsky)
        E[i,j] = SquaBand([kx,ky];t=t,a=a)
    end

    lsE,dos = DOS(E;n = n)
    #dos = (a/2/pi)^2 * dos * dkx
    return lsE,dos,OccupCurve(dos)
end


function HeatCapacity1(Latt::AbstractLattice,lsβ::Union{Vector,StepRangeLen},μ::Number;t=1,a=1)
    L = size(Latt)
    lsk = (1:L)*(pi/(L+1))
    ce = zeros(length(lsβ))
    for (βi,β) in enumerate(lsβ),k in lsk
        ϵk = SquaBand(k;t=t,a=a)
        ce[βi] += (β^2/2/L)*(ϵk-μ)*ϵk/(1+cosh(β*(ϵk-μ)))
    end

    return ce
end

function FreeEnergy1(Latt::AbstractLattice,lsβ::Union{Vector,StepRangeLen},μ::Number;t=1,a=1)
    L = size(Latt)
    lsk = (1:L)*(pi/(L+1))
    f = zeros(length(lsβ))
    for (βi,β) in enumerate(lsβ),k in lsk
        ϵk = SquaBand(k;t=-t,a=a)
        f[βi] += -(1/β)*log(1+exp(-β*(ϵk))) / L
    end

    return f + ParticleNumber1(Latt,lsβ,μ)
end

function InternalEnergy1(Latt::AbstractLattice,lsβ::Union{Vector,StepRangeLen},μ::Number;t=1,a=1)
    L = size(Latt)
    lsk = (1:L)*(pi/(L+1))
    u = zeros(length(lsβ))
    for k in lsk
        ϵk = SquaBand(k;t=-t,a=a)
        for (βi,β) in enumerate(lsβ)
            u[βi] += ϵk/(1+exp(β*(ϵk-μ))) / L
        end
    end

    return u
end

function ParticleNumber1(Latt::AbstractLattice,lsβ::Union{Vector,StepRangeLen},μ::Number;t=1,a=1)
    L = size(Latt)
    lsk = (1:L)*(pi/(L+1))
    n = zeros(length(lsβ))
    for k in lsk
        ϵk = SquaBand(k;t=t,a=a)
        for (βi,β) in enumerate(lsβ)
            n[βi] += 1/(1+exp(β*(ϵk-μ))) / L
        end
    end

    return n
end

