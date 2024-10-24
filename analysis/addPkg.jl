using Pkg

Pkg.activate(".")

Pkg.add("TensorKit")
Pkg.add("JLD2")
Pkg.add("CairoMakie")
Pkg.resolve()
Pkg.gc()

