# iMPS.jl

Personal code for finite MPS simulations, **which is still under great development**.

## Features

### Versatility

Now this code contains:

1. Density matrix renormalization group ([DMRG](https://en.wikipedia.org/wiki/Density_matrix_renormalization_group)) for calculating groundstates.
2. Time-dependent variational principle ([TDVP](https://link.aps.org/doi/10.1103/PhysRevB.94.165116)) for real and imaginary time evolution, which serves for dynamics (i.e., spectrum function and structure factor) and finite temperature simulations (i.e., tangent space tensor renormalization group, [tanTRG](https://link.aps.org/doi/10.1103/PhysRevLett.130.226502)), respectively.
3. Series-expansion thermal tensor network ([SETTN](https://link.aps.org/doi/10.1103/PhysRevB.95.161104)) for finite temperature simulations.

which has been benchmarked and demonstrates a relatively good performance.

### Tutrorials

This code is applied to many models in 1D and 2D finite lattice (square lattice mostly), including:

* Spin: Ising, Heisenberg and XY model.
* Fermion: Free fermion and Hubbard model.

The tutorial is to be added.

### Announcement

Some algorithms in this code is developed by [ CQM2 group](https://www.cqm2itp.com/) which works on purification-based finite-temperature simulations. Relevant packages are listed as follows:

* [FiniteMPS.jl](https://github.com/Qiaoyi-Li/FiniteMPS.jl.git)

## TODO

* Non-abelian symmetries is to be added to accelerate the computation.
* Controlled bond expansion ([CBE](https://doi.org/10.1103/PhysRevLett.130.246402)) is to be implemented to reduce the complexity (of DMRG, TDVP, tanTRG, etc.).
* The code structure is to be optimized to reach a higher performance.

## Acknowledgments

The following packages has been used:

* [TensorKit.jl](https://github.com/Jutho/TensorKit.jl.git) for basic tensor operations.
* [MKL.jl](https://github.com/JuliaLinearAlgebra/MKL.jl.git) for nested multi-threaded BLAS.
