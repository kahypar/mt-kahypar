include("utils.jl")
include("tree_partition/tree_partition.jl")
include("tree_partition/triton_part_refine.jl")
include("overlay.jl")

function tree_distill(embedding::AbstractArray{Float64, 2},
        hgr::__hypergraph__,
        adj_matrix::SparseMatrixCSC,
        hint::AbstractArray{Int},
        hgr_file_in::Union{String, Nothing} = nothing)
    
    hgr_file = isnothing(hgr_file_in) ? write_hypergraph(hgr) : hgr_file_in
    fixed_vertices = __pindex__(Int[], Int[])
    total_weight = sum(hgr.vwts)
    ub_factor = floor(Int, config_e * 50.)
    max_part_weight = Int(max(floor(Int, total_weight * 0.5 * (1.0 + config_e)), ceil(total_weight / 2.)))
    # tree partition
    t = @elapsed begin
        tree_partitions = tree_partition(adj_matrix, 
                                    embedding, 
                                    hgr, 
                                    fixed_vertices, 
                                    ub_factor,
                                    [total_weight - max_part_weight, max_part_weight])
    end
    inform("tree partioning took $(t)s")
    
    partitions = []
    cutsizes = Int[]
    cut_dictionary = Dict{Int, Int}()

    # refine tree partitions
    for i in 1 : length(tree_partitions)
        push!(partitions, triton_part_refine(hgr_file, tree_partitions[i][1], config_k, ub_factor, config_seed, i))
        (cutsize, balance) = golden_evaluator_glob(hgr, config_k, partitions[i])
        push!(cutsizes, cutsize)
        if haskey(cut_dictionary, cutsize) == false
            push!(cut_dictionary, cutsize => i)
        end
        inform("[specpart] Refined partition $i with cutsize $cutsize $balance")
    end

    # select best partitions
    unique_partitions = []
    unique_cutsizes = Int[]
    unique_keys = collect(keys(cut_dictionary))
    for i in 1:length(unique_keys)
        key = unique_keys[i]
        partition_id = cut_dictionary[key]
        push!(unique_cutsizes, key)
        push!(unique_partitions, partitions[partition_id])
    end

    sorted_partition_ids = sortperm(unique_cutsizes)
    best_partitions = [] 
    solns_to_pick = min(config_numSolutions, length(unique_cutsizes))
    for i in 1:solns_to_pick
        push!(best_partitions, unique_partitions[sorted_partition_ids[i]])
        inform("[specpart] partition picked with cutsize $(unique_cutsizes[sorted_partition_ids[i]])")
    end
    push!(best_partitions, hint)

    (solution, cut) = overlay_partitions(best_partitions, hgr, hgr_file, unique_cutsizes[sorted_partition_ids[1]])

    inform("cut from tree distilling $cut")

    if isnothing(hgr_file_in)
        run(`rm -f $hgr_file`)
    end

    return solution, cut
end