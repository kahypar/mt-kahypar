set(UtilSources
        memory_tree.cpp
        )

foreach(partitioner_lib ${PARTITION_LIBRARIES})
  target_sources(${partitioner_lib} PRIVATE ${UtilSources})
endforeach()