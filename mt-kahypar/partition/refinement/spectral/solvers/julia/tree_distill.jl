include("utils.jl")

include("config.jl")

function tree_distill(embedding::AbstractArray{Float64, 2}, hgr::__hypergraph__)
    return length(embedding) < hgr.num_vertices ? embedding[end] : embedding
end