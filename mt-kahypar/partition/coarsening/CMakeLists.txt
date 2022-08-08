set(CoarseningSources
        deterministic_multilevel_coarsener.cpp
        multilevel_uncoarsener.cpp
        nlevel_uncoarsener.cpp
        )

foreach(partitioner_lib ${PARTITION_LIBRARIES})
  target_sources(${partitioner_lib} PRIVATE ${CoarseningSources})
endforeach()
