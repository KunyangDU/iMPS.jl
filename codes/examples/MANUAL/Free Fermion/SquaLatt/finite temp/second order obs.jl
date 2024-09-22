include("../../../../../src/iMPS.jl")
include("../../model.jl")

Lx = 12
Ly = 1

t = 1
μ = 1.0
d = 2

D_MPO = 2^4

Latt = load("examples/MANUAL/Free Fermion/data/finite temp/$(Lx)x$(Ly)/Latt_$(Lx)x$(Ly).jld2")["Latt"]
lsμ = load("examples/MANUAL/Free Fermion/data/finite temp/$(Lx)x$(Ly)/lsμ_$(Lx)x$(Ly)_D_MPO=$(D_MPO).jld2")["lsμ"]
for μ in [1.0]
lsβ1 = load("examples/MANUAL/Free Fermion/data/finite temp/$(Lx)x$(Ly)/lsβ1_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_μ=$(μ).jld2")["lsβ1"]
ni = load("examples/MANUAL/Free Fermion/data/finite temp/$(Lx)x$(Ly)/ni_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_μ=$(μ).jld2")["ni"]
fi = load("examples/MANUAL/Free Fermion/data/finite temp/$(Lx)x$(Ly)/fi_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_μ=$(μ).jld2")["fi"]
ui = load("examples/MANUAL/Free Fermion/data/finite temp/$(Lx)x$(Ly)/ui_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_μ=$(μ).jld2")["ui"]

lsβ2 = centralize(lsβ1)
ce = - lsβ2 .* diff(ui) ./ diff(log.(lsβ1))

@save "examples/MANUAL/Free Fermion/data/finite temp/$(Lx)x$(Ly)/lsβ2_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_μ=$(μ).jld2" lsβ2
@save "examples/MANUAL/Free Fermion/data/finite temp/$(Lx)x$(Ly)/ce_$(Lx)x$(Ly)_D_MPO=$(D_MPO)_μ=$(μ).jld2" ce
end
