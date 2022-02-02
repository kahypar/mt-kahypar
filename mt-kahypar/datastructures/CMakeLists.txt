set(DatastructureSources
        contraction_tree.cpp
        dynamic_hypergraph.cpp
        dynamic_hypergraph_factory.cpp
        dynamic_graph.cpp
        dynamic_graph_factory.cpp
        static_hypergraph_factory.cpp
        static_hypergraph.cpp
        static_graph_factory.cpp
        static_graph.cpp
        graph.cpp
        incident_net_array.cpp
        dynamic_adjacency_array.cpp)

foreach(modtarget IN LISTS TARGETS_WANTING_ALL_SOURCES)
    target_sources(${modtarget} PRIVATE ${DatastructureSources})
endforeach()
