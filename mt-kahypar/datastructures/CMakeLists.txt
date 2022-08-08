set(DatastructureSourcesDefault
      static_hypergraph_factory.cpp
      static_hypergraph.cpp)

set(DatastructureSourcesQuality
      dynamic_hypergraph.cpp
      dynamic_hypergraph_factory.cpp
      incident_net_array.cpp
      contraction_tree.cpp)

set(DatastructureSourcesGraph
      static_graph_factory.cpp
      static_graph.cpp)

set(DatastructureSourcesGraphQuality
      contraction_tree.cpp
      dynamic_graph.cpp
      dynamic_graph_factory.cpp
      dynamic_adjacency_array.cpp)

add_library(MtKaHyParDataStructuresDefault STATIC ${DatastructureSourcesDefault})
add_library(MtKaHyParDataStructuresQuality STATIC ${DatastructureSourcesQuality})
add_library(MtKaHyParDataStructuresGraph STATIC ${DatastructureSourcesGraph})
add_library(MtKaHyParDataStructuresGraphQuality STATIC ${DatastructureSourcesGraphQuality})

list(APPEND DATASTRUCTURE_LIBRARIES
      MtKaHyParDataStructuresDefault
      MtKaHyParDataStructuresQuality
      MtKaHyParDataStructuresGraph
      MtKaHyParDataStructuresGraphQuality)

foreach(datastructure_lib ${DATASTRUCTURE_LIBRARIES})
  set_property(TARGET ${datastructure_lib} PROPERTY CXX_STANDARD 17)
  set_property(TARGET ${datastructure_lib} PROPERTY CXX_STANDARD_REQUIRED ON)
endforeach()

foreach(modtarget IN LISTS TARGETS_WANTING_ALL_SOURCES)
  get_target_property(compile_defs ${modtarget} COMPILE_DEFINITIONS)
  SET(STRONG FALSE)
  SET(GRAPH FALSE)
  SET(GRAPH_TEST FALSE)
  if(${modtarget} MATCHES ".*graph_tests.*")
    SET(GRAPH_TEST TRUE)
  endif()
  foreach( property ${compile_defs} )
    if(${property} MATCHES ".*STRONG.*")
      SET(STRONG TRUE)
    elseif(${property} MATCHES ".*GRAPH.*")
      SET(GRAPH TRUE)
    endif()
  endforeach()

  if ( GRAPH_TEST )
    target_link_libraries(${modtarget} MtKaHyParDataStructuresGraph MtKaHyParDataStructuresGraphQuality)
  elseif ( NOT GRAPH AND NOT STRONG )
    target_link_libraries(${modtarget} MtKaHyParDataStructuresDefault)
  elseif( GRAPH AND NOT STRONG )
    target_link_libraries(${modtarget} MtKaHyParDataStructuresGraph)
  elseif( NOT GRAPH AND STRONG )
    target_link_libraries(${modtarget} MtKaHyParDataStructuresQuality)
  elseif( GRAPH AND STRONG )
    target_link_libraries(${modtarget} MtKaHyParDataStructuresGraphQuality)
  endif()
endforeach()
