set(ProcessMappingSources
        process_graph.cpp
        all_pair_shortest_path.cpp
        steiner_tree.cpp
        process_graph_partitioner.cpp
        )

foreach(modtarget IN LISTS PARTITIONING_SUITE_TARGETS)
    target_sources(${modtarget} PRIVATE ${ProcessMappingSources})
endforeach()