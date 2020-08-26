set(DatastructureSources
        static_hypergraph_factory.cpp
        static_hypergraph.cpp
        graph.cpp)

foreach(modtarget IN LISTS TARGETS_WANTING_ALL_SOURCES)
    target_sources(${modtarget} PRIVATE ${DatastructureSources})
endforeach()
