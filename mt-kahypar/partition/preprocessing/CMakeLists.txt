set(PreprocessingSources
        community_detection/parallel_louvain.cpp
        community_detection/local_moving_modularity.cpp
        ${PROJECT_SOURCE_DIR}/mt-kahypar/datastructures/graph.cpp)

foreach(partitioner_lib ${PARTITION_LIBRARIES})
  target_sources(${partitioner_lib} PRIVATE ${PreprocessingSources})
endforeach()