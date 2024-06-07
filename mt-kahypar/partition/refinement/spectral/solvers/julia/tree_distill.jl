include("utils.jl")
include("config.jl")
include("tree_partition/tree_partition.jl")
include("tree_partition/triton_part_refine.jl")

function tree_distill(embedding::AbstractArray{Float64, 2}, hgr::__hypergraph__, adj_matrix::SparseMatrixCSC, hint::AbstractArray{Int})
    fixed_vertices = __pindex__(Int[], Int[])
    total_weight = sum(hgr.vwts)
    max_part_weight = round(Int, (1 + config_e) * (total_weight * 0.5))
    ub_factor = convert(Int, config_e * 100)
    # tree partition
    tree_partitions = tree_partition(adj_matrix, 
                                    embedding, 
                                    hgr, 
                                    fixed_vertices, 
                                    ub_factor,
                                    [total_weight - max_part_weight, max_part_weight])
    
    partitions = []
    cutsizes = Int[]
    cut_dictionary = Dict{Int, Int}()

    hgr_file = "$config_tmpDir/$(hgr.num_vertices)_$(hgr.num_hyperedges)_$(time()).hgr"
    write_hypergraph(hgr, hgr_file)

    # refine tree partitions
    for i in 1 : length(tree_partitions)
        push!(partitions, triton_part_refine(hgr_file, tree_partitions[i][1], config_k, ub_factor, config_seed, i))
        (cutsize, balance) = golden_evaluator(hgr, config_k, partitions[i])
        push!(cutsizes, cutsize)
        if haskey(cut_dictionary, cutsize) == false
            push!(cut_dictionary, cutsize => i)
        end
        @info "[specpart] Refined partition $i with cutsize $cutsize $balance" 
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
        @info "[specpart] partition picked with cutsize $(unique_cutsizes[sorted_partition_ids[i]])" 
    end
    push!(best_partitions, hint)

    candidate = best_partitions[1]

    run(`rm $hgr_file`)

    return candidate
end