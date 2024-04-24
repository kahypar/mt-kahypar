module CombinatorialMultigrid
    using SparseArrays
    using LinearAlgebra
    using LDLFactorizations
    #using Laplacians

    include("cmg_alg.jl")
    export cmg_preconditioner_lap, lPreconditioner
end