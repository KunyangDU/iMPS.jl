#= using TensorKit,JLD2,MKL,FiniteLattices

include("Environment.jl")
include("Initialize.jl")
include("PushLeft.jl")
include("PushRight.jl")

include("LocalOperator.jl")
include("Node.jl")
include("addIntr1.jl")
include("addIntr2.jl")
include("Automata.jl")
include("Fermion.jl")
include("Spin.jl")

include("MPO.jl")
include("canonicalize.jl")
include("Operation.jl")
include("MPS.jl")
include("Hamiltonian.jl")

include("CalculateObs.jl")
include("methods.jl")

include("Action.jl")
include("Merge.jl")
include("Move.jl")
include("SVD.jl")
include("Variation.jl")
include("Algebra.jl")
include("Geometry.jl")
include("Show.jl")
include("Tools.jl")

include("DMRG.jl")
include("Lanczos.jl")
include("TDVP.jl")
include("SETTN.jl")
include("tanTRG.jl")
 =#

module iMPS
using TensorKit,MKL,FiniteLattices

export PushLeft,PushRight
include("Environment/Environment.jl")
include("Environment/Initialize.jl")
include("Environment/PushLeft.jl")
include("Environment/PushRight.jl")

export AbstractLocalOperator,IdentityOperator,LocalOperator
export InteractionTreeNode,InteractionTree
export addIntr1!,addIntr2!,AutomataMPO
include("IntrTree/LocalOperator.jl")
include("IntrTree/Node.jl")
include("IntrTree/addIntr1.jl")
include("IntrTree/addIntr2.jl")
include("IntrTree/Automata.jl")

export Spin2,Spin2Fermion,SpinlessFermion
include("LocalSpace/Fermion.jl")
include("LocalSpace/Spin.jl")

export RandMPS
include("MatProd/MPO/MPO.jl")
include("MatProd/MPO/canonicalize.jl")
include("MatProd/MPO/Operation.jl")
include("MatProd/MPS/MPS.jl")
include("MatProd/Hamiltonian.jl")

export calObs
include("Observables/CalculateObs.jl")
include("Observables/methods.jl")

export showQuantSweep,showBlockMPO,showdomain,FindMaxDist
include("Operation/Action.jl")
include("Operation/Merge.jl")
include("Operation/Move.jl")
include("Operation/SVD.jl")
include("Operation/Variation.jl")
include("Tools/Algebra.jl")
include("Tools/Geometry.jl")
include("Tools/Show.jl")
include("Tools/Tools.jl")

export sweepDMRG2
include("Algorithm/DMRG.jl")
include("Algorithm/Lanczos.jl")
include("Algorithm/TDVP.jl")
include("Algorithm/SETTN.jl")
include("Algorithm/tanTRG.jl")
end

#= 
MPO data matrix should be the hermitian conjugate of the 
matrix of observables, i.e.
M_{data} = M_{Obs}⁺
under this convention, the evolution operator is exp(i)
=#

#= 
Suppose you have been familiar with TensorKit.jl
Indexing formalism of TensorMap
[] -> codomain, () -> domain

I. MPS
    1. Right Orthogonal:
           ___
    (2) → |   | → [1]
           ‾‾‾
            ↑
           (3)
    i.e. TensorMap{ [1] , (2) ⊗ (3) }
           [2] 
            ↑
           ___
    [1] ← |   | ← (3)
           ‾‾‾

    i.e. TensorMap{ [1] ⊗ [2], (3)  }
    -------------------
    2. Left Orthogonal:
           ___
    [1] ← |   | ← (3)
           ‾‾‾
            ↑
           (2)
    i.e. TensorMap{ [1] , (2) ⊗ (3) }
           [1]
            ↑
           ___
    (3) → |   | → [2]
           ‾‾‾
    i.e. TensorMap{ [1] ⊗ [2]  , (3) }
    -------------------
    3. Center:
           ___
    (1) → |   | ← (3)
           ‾‾‾
            ↑
           (2)
    i.e. TensorMap{ (1) ⊗ (2) ⊗ (3)}
           [2]
            ↑
           ___
    [1] ← |   | → [3]
           ‾‾‾
    i.e. TensorMap{ [1] ⊗ [2] ⊗ [3]}
    -------------------
    4. SVD conjunction
    -Right
           ___
    (2) → | → | ← (1)
           ‾‾‾
    i.e. TensorMap{ (1) ⊗ (2) }

    -Left
           ___
    (1) → | ← | ← (2)
           ‾‾‾
    i.e. TensorMap{ (1) ⊗ (2) }
    -------------------
    5. 2-site Center MPS
           _______
    (1) → |       | ← (4)
           ‾‾‾‾‾‾‾
            ↑   ↑
           (2) (3)
    i.e. TensorMap{ (1) ⊗ (2) ⊗ (3) ⊗ (4) }
           [2] [3]
            ↑   ↑
           _______
    [1] ← |       | → [4]
           ‾‾‾‾‾‾‾
    i.e. TensorMap{ [1] ⊗ [2] ⊗ [3] ⊗ [4] }

II.MPO
    1.1 Local MPO (Right Orthogonal)
           [2]
            ↑
           ___
    (3) → |   | → [1]
           ‾‾‾
            ↑
           (4)
    i.e. TensorMap{ [1] ⊗ [2], (3) ⊗ (4) }

           (4)
            ↓
           ___
    [1] ← |   | ← (3)
           ‾‾‾
            ↓
           [2]
    i.e. TensorMap{ [1] ⊗ [2], (3) ⊗ (4) }

    1.2 Local MPO (Left Orthogonal)
           [1]
            ↑
           ___
    [2] ← |   | ← (4)
           ‾‾‾
            ↑
           (3)
    i.e. TensorMap{ [1] ⊗ [2], (3) ⊗ (4) }

           (3)
            ↓
           ___
    (4) → |   | → [2]
           ‾‾‾
            ↓
           [1]
    i.e. TensorMap{ [1] ⊗ [2], (3) ⊗ (4) }

    1.3 Local MPO (Center Orthogonal)
           [1]
            ↑
           ___
    (2) → |   | ← (4)
           ‾‾‾
            ↑
           (3)
    i.e. TensorMap{ [1], (2) ⊗ (3) ⊗ (4) }

           (4)
            ↓
           ___
    [1] ← |   | → [3]
           ‾‾‾
            ↓
           [2]
    i.e. TensorMap{ [1] ⊗ [2] ⊗ [3], (4) }

    1.4 Local 2 site MPO (Center Orthogonal)
           [2] [1]
            ↑   ↑
           _______
    (3) → |       | ← (6)
           ‾‾‾‾‾‾‾
            ↑   ↑
           (4) (5)
    i.e. TensorMap{ [1] ⊗ [2], (3) ⊗ (4) ⊗ (5) ⊗ (6) }

           (4)
            ↓
           ___
    [1] ← |   | → [3]
           ‾‾‾
            ↓
           [2]
    i.e. TensorMap{ [1] ⊗ [2] ⊗ [3], (4) }

    2. 0-site Effective Hamiltonian
         [1] [2]
          ↑   ↑
         _______
        |       |
         ‾‾‾‾‾‾‾
          ↑   ↑
         (3) (4)
    i.e. TensorMap{ [1] ⊗ [2], (3) ⊗ (4) }

    3. 1-site Effective Hamiltonian
         [1] [2] [3]
          ↑   ↑   ↑
         ___________
        |           |
         ‾‾‾‾‾‾‾‾‾‾‾
          ↑   ↑   ↑
         (4) (5) (6)
    i.e. TensorMap{ [1] ⊗ [2] ⊗ [3], (4) ⊗ (5) ⊗ (6) }

    4. 2-site Effective Hamiltonian
         [1] [2] [3] [4]
          ↑   ↑   ↑   ↑
         _______________
        |               |
         ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
          ↑   ↑   ↑   ↑
         (5) (6) (7) (8)
    i.e. TensorMap{ [1] ⊗ [2] ⊗ [3] ⊗ [4], (5) ⊗ (6) ⊗ (7) ⊗ (8) }

III.Environment
    1. Right Environment
           ___
    [1] ← |   |
    (2) → |   |
    (3) → |   |
           ‾‾‾
    i.e. TensorMap{ [1] , (2) ⊗ (3) }

    *canonicalized
           ___
    [1] ← |   |
    [2] ← |   |
    (3) → |   |
           ‾‾‾
    i.e. TensorMap{ [1] ⊗ [2], (3) }

           ___
    [1] ← |   |
    (2) → |   |
           ‾‾‾
    i.e. TensorMap{ [1] , (2) }

    2. Left Environment
           ___
          |   | → [1]
          |   | → [2]
          |   | ← (3)
           ‾‾‾
    i.e. TensorMap{ [1] ⊗ [2], (3)  }
           ___
          |   | → [1]
          |   | ← (2)
           ‾‾‾
    i.e. TensorMap{ [1] , (2)  }

 =#



