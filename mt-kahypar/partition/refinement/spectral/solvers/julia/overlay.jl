include("overlay/contraction.jl")

function overlay_partitions(partitions, hgr)
    hgraph_contracted, clusters = overlay(best_partitions, hgraph)
    refined_partition = optimal_partitioner(hmetis_path, ilp_path, hgraph_contracted, num_parts, ub_factor)
    cutsize = golden_evaluator(hgraph_contracted, num_parts, refined_partition)
    partition_projected = zeros(Int, hgraph.num_vertices)
    for i in 1:length(clusters)
        cid = clusters[i]
        partition_projected[i] = refined_partition[cid]
    end
    specpart_partition_name = processed_hg_name * ".specpart" * ".part.2"
    write_partition(partition_projected, specpart_partition_name)
    triton_part_refine(triton_part_refiner_path, processed_hg_name, specpart_partition_name, num_parts, ub_factor, seed, 0)
    partition_projected = read_hint_file(specpart_partition_name)
    cutsize = golden_evaluator(hgraph, num_parts, partition_projected)
    
    return partition_projected, cutsize
end