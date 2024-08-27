# adapted from [K-SpecPart](https://github.com/TILOS-AI-Institute/HypergraphPartitioning/tree/main/K_SpecPart) under [BSD license](https://github.com/TILOS-AI-Institute/HypergraphPartitioning/blob/main/LICENSE)

module rmq

using LightGraphs
using SimpleWeightedGraphs

include("euler_tour.jl")
include("node_levels.jl")
include("rmq_solve.jl")
include("queries.jl")
include("lca2rmq.jl")

export lca2rmq

end