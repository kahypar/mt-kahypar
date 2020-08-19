set(IOSources csv_output.cpp hypergraph_io.cpp sql_plottools_serializer.cpp)

foreach(modtarget IN LISTS TARGETS_WANTING_ALL_SOURCES)
    target_sources(${modtarget} PRIVATE ${IOSources})
endforeach()
