set(InitialPartitioningSources
        flat/pool_initial_partitioner.cpp
        flat/random_initial_partitioner.cpp
        flat/bfs_initial_partitioner.cpp
        flat/label_propagation_initial_partitioner.cpp
        recursive_bisection_initial_partitioner.cpp
        recursive_initial_partitioner.cpp
        )

foreach(modtarget IN LISTS TARGETS_WANTING_ALL_SOURCES)
    target_sources(${modtarget} PRIVATE ${InitialPartitioningSources})
endforeach()