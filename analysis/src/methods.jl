

function centralize(data::Union{Vector,OrdinalRange,StepRangeLen})
    return [(data[i] + data[i+1]) / 2 for i in 1:(length(data) - 1)]
end

function DOS(E::Array;n::Int64=size(E)[1])
    secs = range(extrema(E)...,n)
    dos = zeros(Int64,n-1)
    lsE = centralize(secs)

    for i in eachindex(secs)[1:end-1]
        dos[i] = count(x -> secs[i] <= x <= secs[i+1], E)
    end

    return lsE,dos / n^2
end

function OccupCurve(dos::Vector)
    return cumsum(dos)
end

