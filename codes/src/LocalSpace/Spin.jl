module Spin2
using TensorKit
phys = (ℂ^2)'

const Spinσ0 = TensorMap([1 0;0 1],phys,phys)
const Spinσx = TensorMap([0 1;1 0],phys,phys)
const Spinσy = TensorMap([0 -1im;1im 0],phys,phys)
const Spinσz = TensorMap([1 0; 0 -1],phys,phys)

const I = let 
    Spinσ0
end

const Sx = let 
    Spinσx / 2
end

const Sy = let 
    Spinσy / 2
end

const Sz = let 
    Spinσz / 2
end

const SxSx = let 
    (Spinσx,Spinσx) ./ 2
end

const SySy = let 
    (Spinσy,Spinσy) ./ 2
end

const SzSz = let 
    (Spinσz,Spinσz) ./ 2
end

const S₊S₋ = let 
    (Spinσx + 1im * Spinσx) / 4, (Spinσx - 1im * Spinσx) / 4
end

const S₋S₊ = let 
    (Spinσx - 1im * Spinσx) / 4, (Spinσx + 1im * Spinσx) / 4
end

end