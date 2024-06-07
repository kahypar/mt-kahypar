using LinearAlgebra
using LinearMaps
using IterativeSolvers
using SparseArrays
using Random

include("config.jl")
include("utils.jl")
include("graph/hypergraph.jl")
include("graph/matrices.jl")
include("lobpcg.jl")
include("tree_distill.jl")

function main(input, method)
    inform("julia launched")
    n = convert(Int, input[1][1])
    inform_dbg(n, false, () -> "transmitted (hyper)graph data: " * string(convert(AbstractArray{Int64}, input[1])))

    try
        hgr = import_hypergraph(input[1])
        hint_partition = convert(AbstractArray{Int64}, input[2])
        deflation_space = reshape(convert(AbstractArray{Float64}, input[3]), n, convert(Int, length(input[3]) / n))

        inform_dbg(n, true, "building adjaciency matrix...")
        adj_matrix = adjacency_matrix(hgr)
        method(hgr, hint_partition, deflation_space, adj_matrix)
    catch e
        inform("failed due to " * sprint(showerror, e))
        @print_backtrace
        
        return ones(Float64, n)
    end
end

function main_lobpcg(hgr_data::AbstractArray, hint::AbstractArray, deflation_evecs::AbstractArray)
    main((hgr_data, hint, deflation_evecs), (g, h, e, a) -> solve_lobpcg(g, h, e, 1, a))
end

function main_lobpcg_tree_distill(hgr_data::AbstractArray, hint::AbstractArray, deflation_evecs::AbstractArray)
    main((hgr_data, hint, deflation_evecs), function(g::__hypergraph__, h::AbstractArray{Int}, e::AbstractArray{Float64, 2}, a::SparseMatrixCSC)
        return tree_distill(solve_lobpcg(g, h, e, config_numEvecs, a), g, a, h) 
    end)
end

main_auto = getfield(Main, Symbol(config_main))
