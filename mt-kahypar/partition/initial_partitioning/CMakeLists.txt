set(InitialPartitioningSources
        flat/pool_initial_partitioner.cpp
        flat/random_initial_partitioner.cpp
        flat/bfs_initial_partitioner.cpp
        flat/label_propagation_initial_partitioner.cpp
        recursive_bipartitioning_initial_partitioner.cpp
        deep_initial_partitioner.cpp
        )

foreach(partitioner_lib ${PARTITION_LIBRARIES})
  target_sources(${partitioner_lib} PRIVATE ${InitialPartitioningSources})
endforeach()