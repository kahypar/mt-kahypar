# adapted from [K-SpecPart](https://github.com/TILOS-AI-Institute/HypergraphPartitioning/tree/main/K_SpecPart) under [BSD license](https://github.com/TILOS-AI-Institute/HypergraphPartitioning/blob/main/LICENSE)

using SimpleWeightedGraphs

include("../graph/hypergraph.jl")

struct __pindex__
    p1::Vector{Int}
    p2::Vector{Int}
end

struct __least_common_ancestor__
    rmq_sparse_table::Matrix{Int}
    euler_level::Vector{Int}
    child::Vector{Int}
    parents::Vector{Int}
    euler_tour::Vector{Int}
    level_vec::Vector{Int}
    fts::Vector{Int}
    ifts::Vector{Int}
end

struct __cut_profile__
    vtx_cuts::AbstractArray{Int}
    edge_cuts::Vector{Int}
    edge_diff::Vector{Int}
    pred::Vector{Int}
    edge_terminators::Vector{Int}
    p::__pindex__
    forced_type::Vector{Int}
    forced_0::Vector{Int}
    forced_1::Vector{Int}
    forced_01::Vector{Int}
    FB0::Vector{Int}
    FB1::Vector{Int}
    edge_cuts_0::Vector{Int}
    edge_cuts_1::Vector{Int}
end

struct __best_partition__
    total_cost::Float64
    area_cost::Float64
    cut_cost::Float64
    cutsize::Int
    partition::Vector{Int}
    cut_point::Int
    area::Vector{Int}
end

struct __recursive_parts__
    hypergraph::__hypergraph__
    T::SimpleWeightedGraphs.SimpleGraph
    distilled_cuts::__cut_profile__
    capacities::Vector{Int}
    cluster_labels::Vector{Int}
end

mutable struct __tree_cuts__
    nforced0::Int
    nforced1::Int
    nforced01::Int
    total_vwts::Int
    exc0::Vector{Int}
    exc1::Vector{Int}
    area::Vector{Int}
    cut_cost0::Vector{Float64}
    cut_cost1::Vector{Float64}
    ratio_cost::Vector{Float64}
    cut_cost::Vector{Float64}
    area_cost::Vector{Float64}
    total_cost::Vector{Float64}
    polarity::Vector{Int}
    status_flag::Vector{Int}
    area_util0::Vector{Int}
    area_util1::Vector{Int}
    pred::Vector{Int}
    hyperedges_flag::Vector{Int}
end