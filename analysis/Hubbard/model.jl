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



