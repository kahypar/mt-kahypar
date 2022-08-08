set(IOSources
        csv_output.cpp
        hypergraph_io.cpp
        sql_plottools_serializer.cpp
        partitioning_output.cpp
        command_line_options.cpp)

foreach(partitioner_lib ${PARTITION_LIBRARIES})
  target_sources(${partitioner_lib} PRIVATE ${IOSources})
endforeach()
