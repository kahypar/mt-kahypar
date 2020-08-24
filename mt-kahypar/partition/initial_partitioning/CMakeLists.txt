set(InitialPartitioningSources
        flat/pool_initial_partitioner.cpp
        recursive_bisection_initial_partitioner.cpp
        )

foreach(modtarget IN LISTS TARGETS_WANTING_ALL_SOURCES)
    target_sources(${modtarget} PRIVATE ${InitialPartitioningSources})
endforeach()