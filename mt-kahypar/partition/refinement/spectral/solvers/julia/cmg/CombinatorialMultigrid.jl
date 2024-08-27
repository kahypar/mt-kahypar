# adapted from [K-SpecPart](https://github.com/TILOS-AI-Institute/HypergraphPartitioning/tree/main/K_SpecPart) under [BSD license](https://github.com/TILOS-AI-Institute/HypergraphPartitioning/blob/main/LICENSE)

module CombinatorialMultigrid
    using SparseArrays
    using LinearAlgebra
    using LDLFactorizations
    #using Laplacians

    include("cmg_alg.jl")
    export cmg_preconditioner_lap, lPreconditioner
end