using Pkg

Pkg.activate(".")

Pkg.add("TensorKit")
Pkg.add("JLD2")
Pkg.add("MKL")
Pkg.add("BenchmarkTools")

Pkg.resolve()
Pkg.gc()

