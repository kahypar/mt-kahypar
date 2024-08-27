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


# generic main method, not to be called directly
function main(input, method, rescue=true)
    inform("julia launched")
    n = convert(Int, input[1][1])
    inform_dbg(n, false, () -> "transmitted (hyper)graph data: " * string(convert(AbstractArray{Int64}, input[1])))

    try
        hgr = import_hypergraph(input[1])
        (hgr_processed, original_indices, new_indices, unused_indices) = isolate_islands(hgr)
        n_processed = hgr_processed.num_vertices
        inform("removed islands containing $(n - n_processed) nodes")

        hint_partition = convert(AbstractArray{Int64}, input[2])
        hint_processed = @view hint_partition[original_indices]

        deflation_space = reshape(convert(AbstractArray{Float64}, input[3]), n, convert(Int, length(input[3]) / n))
        deflation_space_processed = deflation_space[original_indices, :]

        adj_matrix = adjacency_matrix(hgr_processed)

        global global_best = nothing
        golden_evaluator_glob(hgr_processed, config_k, hint_processed)

        # apply chosen algo variant
        result_processed = method(hgr_processed, hint_processed, deflation_space_processed, adj_matrix)

        result = zeros(Float64, n)
        result[original_indices] = result_processed
        result[unused_indices] = hint_partition[unused_indices]

        return result
    catch e
        inform("failed due to " * sprint(showerror, e))
        @print_backtrace
        
        if rescue
            return ones(Float64, n)
        else
            exit(1)
        end
    end
end

# only solve gevp with lobpcg and return fiedler vector
function main_lobpcg(hgr_data::AbstractArray, hint::AbstractArray, deflation_evecs::AbstractArray, rescue=true)
    return main((hgr_data, hint, deflation_evecs, rescue), (g, h, e, a) -> solve_lobpcg(g, h, e, 1, a))
end

# solve gevp and then perform tree distilling
function main_lobpcg_tree_distill(hgr_data::AbstractArray, hint::AbstractArray, deflation_evecs::AbstractArray)
    return main((hgr_data, hint, deflation_evecs), (g, h, e, a) -> tree_distill(solve_lobpcg(g, h, e, config_numEvecs, a), g, a, h)[1])
end

# apply whole algorithm adapted from kspecpart
function main_kspecpart(hgr_data::AbstractArray, hint::AbstractArray, deflation_evecs::AbstractArray)
    return main((hgr_data, hint, deflation_evecs), function (g, h, e, a)
        candidates = []
        hgr_file = write_hypergraph(g)
        hint = h
        for i in 1 : config_numCandidates
            push!(candidates, tree_distill(solve_lobpcg(g, hint, e, config_numEvecs, a), g, a, hint, hgr_file))
            hint = candidates[end][1]
        end
        sort!(candidates, by = x->x[2])
        overlay_partitions(map(c->c[1], candidates), g, hgr_file, candidates[1][2])
        inform("best cut from specpart $(global_best[1])")
        return global_best[2]
    end)
end

main_lobpcg_unsafe(g, h, d) = main_lobpcg(g, h, d, false)

main_nothing = (d, h, e) -> ones(Float64, d[1])

# retrieve main variant from config file
main_auto = getfield(Main, Symbol(config_main))

global_best = nothing
