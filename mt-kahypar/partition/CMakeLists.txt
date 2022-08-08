add_subdirectory(preprocessing/)
add_subdirectory(coarsening)
add_subdirectory(refinement)
add_subdirectory(initial_partitioning)
add_subdirectory(registries)

set(PartitionSources
        partitioner.cpp
        multilevel.cpp
        context.cpp
        context_enum_classes.cpp
        metrics.cpp
        recursive_bipartitioning.cpp
        deep_multilevel.cpp
        )

foreach(partitioner_lib ${PARTITION_LIBRARIES})
  target_sources(${partitioner_lib} PRIVATE ${PartitionSources})
endforeach()