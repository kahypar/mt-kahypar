set(MultilevelDatastructureSources
      static_hypergraph_factory.cpp
      static_hypergraph.cpp
      graph.cpp)

set(MultilevelGraphDatastructureSources
      static_graph_factory.cpp
      static_graph.cpp
      graph.cpp)

set(NLevelDatastructureSources
      contraction_tree.cpp
      dynamic_hypergraph.cpp
      dynamic_hypergraph_factory.cpp
      graph.cpp
      incident_net_array.cpp)

set(NLevelGraphDatastructureSources
      contraction_tree.cpp
      dynamic_graph.cpp
      dynamic_graph_factory.cpp
      graph.cpp
      dynamic_adjacency_array.cpp)

foreach(modtarget IN LISTS TARGETS_WANTING_ALL_SOURCES)
    get_target_property(compile_defs ${modtarget} COMPILE_DEFINITIONS)
    SET(QUALITY FALSE)
    SET(GRAPH FALSE)
    SET(GRAPH_TEST FALSE)
    if(${modtarget} MATCHES ".*graph_tests.*")
      SET(GRAPH_TEST TRUE)
    endif()
    foreach( property ${compile_defs} )
      if(${property} MATCHES ".*QUALITY.*")
        SET(QUALITY TRUE)
      elseif(${property} MATCHES ".*GRAPH.*")
        SET(GRAPH TRUE)
      endif()
    endforeach()

    if ( GRAPH_TEST )
      target_sources(${modtarget} PRIVATE ${MultilevelGraphDatastructureSources} ${NLevelGraphDatastructureSources})
    elseif ( NOT GRAPH AND NOT QUALITY )
      target_sources(${modtarget} PRIVATE ${MultilevelDatastructureSources})
    elseif( GRAPH AND NOT QUALITY )
      target_sources(${modtarget} PRIVATE ${MultilevelGraphDatastructureSources})
    elseif( NOT GRAPH AND QUALITY )
      target_sources(${modtarget} PRIVATE ${NLevelDatastructureSources})
    elseif( GRAPH AND QUALITY )
      target_sources(${modtarget} PRIVATE ${NLevelGraphDatastructureSources})
    endif()
endforeach()
