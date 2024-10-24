include("../../../src/iMPS.jl")
include("../model.jl")

Lx = 8
Ly = 1

t = 1
d = 2

D_MPO = 20

Latt = YCSqua(Lx,Ly)
#lsμ = load("examples/Free Fermion/data/finite temp/$(Lx)x$(Ly)/lsμ_$(Lx)x$(Ly)_D_MPO=$(D_MPO).jld2")["lsμ"]
for μ in [2.0]
lsβ1 = load("examples/Spinless Fermion/data/finite temp/lsβ1_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_t=$(t)_μ=$(μ).jld2")["lsβ1"]
ni = load("examples/Spinless Fermion/data/finite temp/ni_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_t=$(t)_μ=$(μ).jld2")["ni"]
fi = load("examples/Spinless Fermion/data/finite temp/fi_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_t=$(t)_μ=$(μ).jld2")["fi"]
ui = load("examples/Spinless Fermion/data/finite temp/ui_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_t=$(t)_μ=$(μ).jld2")["ui"]

lsβ2 = centralize(lsβ1)
ce = - lsβ2 .* diff(ui) ./ diff(log.(lsβ1))

@save "examples/Spinless Fermion/data/finite temp/lsβ2_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_t=$(t)_μ=$(μ).jld2" lsβ2
@save "examples/Spinless Fermion/data/finite temp/ce_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_t=$(t)_μ=$(μ).jld2" ce
end
