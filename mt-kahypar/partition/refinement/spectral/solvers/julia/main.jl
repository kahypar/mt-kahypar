using LinearAlgebra
using LinearMaps
using IterativeSolvers
using SparseArrays
using Random

include("config.jl")
include("utils.jl")
include("graph/hypergraph.jl")
include("graph/matrices.jl")
include("graph/isolate_islands.jl")
include("lobpcg.jl")
include("tree_distill.jl")

function main(input, method)
    inform("julia launched")
    n = convert(Int, input[1][1])
    inform_dbg(n, false, () -> "transmitted (hyper)graph data: " * string(convert(AbstractArray{Int64}, input[1])))

    try
        hgr = import_hypergraph(input[1])
        (hgr_processed, original_indices, new_indices, unused_indices) = isolate_islands(hgr)
        n_processed = hgr_processed.num_vertices

        hint_partition = convert(AbstractArray{Int64}, input[2])
        hint_processed = @view hint_partition[original_indices]

        deflation_space = reshape(convert(AbstractArray{Float64}, input[3]), n, convert(Int, length(input[3]) / n))
        deflation_space_processed = deflation_space[original_indices, :]

        adj_matrix = adjacency_matrix(hgr_processed)

        result_processed = method(hgr_processed, hint_processed, deflation_space_processed, adj_matrix)

        result = zeros(Float64, n)
        result[original_indices] = result_processed
        result[unused_indices] = hint_partition[unused_indices]

        return result
    catch e
        inform("failed due to " * sprint(showerror, e))
        @print_backtrace
        
        return ones(Float64, n)
    end
end

function main_lobpcg(hgr_data::AbstractArray, hint::AbstractArray, deflation_evecs::AbstractArray)
    return main((hgr_data, hint, deflation_evecs), (g, h, e, a) -> solve_lobpcg(g, h, e, 1, a))
end

function main_lobpcg_tree_distill(hgr_data::AbstractArray, hint::AbstractArray, deflation_evecs::AbstractArray)
    return main((hgr_data, hint, deflation_evecs), (g, h, e, a) -> tree_distill(solve_lobpcg(g, h, e, config_numEvecs, a), g, a, h))
end

function main_kspecpart(hgr_data::AbstractArray, hint::AbstractArray, deflation_evecs::AbstractArray)
    return main((hgr_data, hint, deflation_evecs), function (g, h, e, a)
        candidates = []
        hgr_file = write_hypergraph(g)
        for i in 1 : config_numCandidates
            push!(candidates, tree_distill(solve_lobpcg(g, h, e, config_numEvecs, a), g, a, h, hgr_file))
        end
        overlay_partitions(candidates, g, hgr_file)
    end)
end

main_auto = getfield(Main, Symbol(config_main))
