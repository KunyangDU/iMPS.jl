function GroupV(k::Number,J::Number,h::Number)
    return J*h*sin(k)/sqrt(4*h^2 + J^2 - 4*J*h*cos(k))
end

function MaxGroupV(J::Number,h::Number)
    return GroupV(pi/2,J,h)
end


function IsingBand(k::Number,J::Number,h::Number;a::Number=1)
    return 2*sqrt(J^2 + h^2 + 2*J*h*cos(k*a))
end





