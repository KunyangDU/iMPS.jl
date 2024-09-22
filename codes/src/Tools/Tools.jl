
function irange(start::Int64,stop::Int64)
    if start > stop
        return start:-1:stop 
    else
        return start:stop
    end
end

function FindPair(pairs::Vector,site::Number,target::Number)
    return filter(x -> x[site] == target, pairs)
end

function ApproxReal(Qi::Number;tol::Float64=1e-5)
    imag(Qi) <= tol && return real(Qi)
    @error "not real"
end


function centralize(data::Union{Vector,OrdinalRange,StepRangeLen})
    return [(data[i] + data[i+1]) / 2 for i in 1:(length(data) - 1)]
end





