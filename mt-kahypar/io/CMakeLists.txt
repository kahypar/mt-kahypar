set(IOSources
        csv_output.cpp
        hypergraph_io.cpp
        hypergraph_factory.cpp
        sql_plottools_serializer.cpp
        partitioning_output.cpp)

if(NOT MT_KAHYPAR_DISABLE_BOOST)
    list(APPEND IOSources command_line_options.cpp)
endif()

foreach(modtarget IN LISTS PARTITIONING_SUITE_TARGETS)
    target_sources(${modtarget} PRIVATE ${IOSources})
endforeach()

set(ToolsIOSources
        hypergraph_io.cpp
        hypergraph_factory.cpp
        partitioning_output.cpp)

foreach(modtarget IN LISTS TOOLS_TARGETS)
    target_sources(${modtarget} PRIVATE ${ToolsIOSources})
endforeach()