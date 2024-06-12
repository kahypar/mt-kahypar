include("overlay/contraction.jl")

function overlay_partitions(partitions, hgr, hgr_file)
    ub_factor = ceil(Int, config_e * 100)
    hgraph_contracted, clusters = overlay(partitions, hgr)
    refined_partition = optimal_partitioner(hgraph_contracted, config_k, ub_factor)
    cutsize = golden_evaluator(hgraph_contracted, config_k, refined_partition)
    partition_projected = zeros(Int, hgr.num_vertices)
    for i in 1:length(clusters)
        cid = clusters[i]
        partition_projected[i] = refined_partition[cid]
    end
    specpart_partition_name = hgr_file * ".specpart" * ".part.2"
    write_partition(partition_projected, specpart_partition_name)
    triton_part_refine(hgr_file, specpart_partition_name, config_k, ub_factor, config_seed, 0)
    partition_projected = read_hint_file(specpart_partition_name)
    cutsize = golden_evaluator(hgr, config_k, partition_projected)
    
    return partition_projected, cutsize
end