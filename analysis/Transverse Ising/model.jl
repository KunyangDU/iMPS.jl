function GroupV(k::Number,J::Number,h::Number)
    return 2*J*h*sin(k)/sqrt(h^2 + J^2 - 2*J*h*cos(k))
end

function MaxGroupV(J::Number,h::Number)
    return GroupV(pi/2,J,h)
end