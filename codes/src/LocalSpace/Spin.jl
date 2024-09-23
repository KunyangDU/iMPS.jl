module Spin2
using TensorKit

const d = 2
const PhySpace = (ℂ^d)'
const Spinσ0 = TensorMap([1 0;0 1],PhySpace,PhySpace)
const Spinσx = TensorMap([0 1;1 0],PhySpace,PhySpace)
const Spinσy = TensorMap([0 -1im;1im 0],PhySpace,PhySpace)
const Spinσz = TensorMap([1 0; 0 -1],PhySpace,PhySpace)

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