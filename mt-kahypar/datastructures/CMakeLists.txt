set(DatastructureSources
        contraction_tree.cpp
        dynamic_hypergraph.cpp
        dynamic_hypergraph_factory.cpp
        static_hypergraph_factory.cpp
        static_hypergraph.cpp
        graph.cpp
        incident_net_array.cpp
        async/uncontraction_group_tree.cpp
        async/depth_priority_queue.cpp
        )

foreach(modtarget IN LISTS TARGETS_WANTING_ALL_SOURCES)
    target_sources(${modtarget} PRIVATE ${DatastructureSources})
endforeach()
