set(IOSources
        csv_output.cpp
        hypergraph_io.cpp
        hypergraph_factory.cpp
        sql_plottools_serializer.cpp
        partitioning_output.cpp
        presets.cpp)

if(NOT KAHYPAR_DISABLE_BOOST)
    list(APPEND IOSources command_line_options.cpp)
endif()

target_sources(MtKaHyPar-Sources INTERFACE ${IOSources})

set(ToolsIOSources
        hypergraph_io.cpp
        hypergraph_factory.cpp
        partitioning_output.cpp)

target_sources(MtKaHyPar-ToolsSources INTERFACE ${ToolsIOSources})
