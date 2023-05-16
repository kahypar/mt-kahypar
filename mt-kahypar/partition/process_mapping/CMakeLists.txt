set(ProcessMappingSources
        process_graph.cpp
        all_pair_shortest_path.cpp
        )

foreach(modtarget IN LISTS PARTITIONING_SUITE_TARGETS)
    target_sources(${modtarget} PRIVATE ${ProcessMappingSources})
endforeach()