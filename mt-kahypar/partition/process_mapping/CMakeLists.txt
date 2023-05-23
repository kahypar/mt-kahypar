set(ProcessMappingSources
        process_graph.cpp
        all_pair_shortest_path.cpp
        steiner_tree.cpp
        dual_bipartitioning.cpp
        greedy_mapping.cpp
        initial_mapping.cpp
        kerninghan_lin.cpp
        )

foreach(modtarget IN LISTS PARTITIONING_SUITE_TARGETS)
    target_sources(${modtarget} PRIVATE ${ProcessMappingSources})
endforeach()