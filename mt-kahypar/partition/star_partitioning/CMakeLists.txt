set(SPSources
        star_partitioning.cpp
        simple_greedy.cpp
        approximate.cpp
        rate_separated.cpp
        detect_low_degree_nodes.cpp
        )

foreach(modtarget IN LISTS TARGETS_WANTING_ALL_SOURCES)
    target_sources(${modtarget} PRIVATE ${SPSources})
endforeach()
