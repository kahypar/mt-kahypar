set(DatastructureSources
      graph.cpp
      fixed_vertex_support.cpp)

set(MultilevelDatastructureSources
      static_hypergraph_factory.cpp
      static_hypergraph.cpp)

set(MultilevelGraphDatastructureSources
      static_graph_factory.cpp
      static_graph.cpp)

set(NLevelDatastructureSources
      contraction_tree.cpp
      dynamic_hypergraph.cpp
      dynamic_hypergraph_factory.cpp
      incident_net_array.cpp)

set(NLevelGraphDatastructureSources
      contraction_tree.cpp
      dynamic_graph.cpp
      dynamic_graph_factory.cpp
      dynamic_adjacency_array.cpp)

target_sources(MtKaHyPar-Sources INTERFACE ${DatastructureSources})
target_sources(MtKaHyPar-Sources INTERFACE ${MultilevelDatastructureSources})
if(KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES)
      target_sources(MtKaHyPar-Sources INTERFACE ${MultilevelGraphDatastructureSources})
else()
      target_sources(MtKaHyPar-LibraryBuildSources PRIVATE ${MultilevelGraphDatastructureSources})
endif()
if(KAHYPAR_ENABLE_HIGHEST_QUALITY_FEATURES)
      target_sources(MtKaHyPar-Sources INTERFACE ${NLevelDatastructureSources})
else()
      target_sources(MtKaHyPar-LibraryBuildSources PRIVATE ${NLevelDatastructureSources})
endif()
if(KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES AND KAHYPAR_ENABLE_HIGHEST_QUALITY_FEATURES)
      target_sources(MtKaHyPar-Sources INTERFACE ${NLevelGraphDatastructureSources})
else()
      target_sources(MtKaHyPar-LibraryBuildSources PRIVATE ${NLevelGraphDatastructureSources})
endif()


set(ToolsDatastructureSources
      static_hypergraph_factory.cpp
      static_hypergraph.cpp
      static_graph_factory.cpp
      static_graph.cpp
      contraction_tree.cpp
      dynamic_hypergraph.cpp
      dynamic_hypergraph_factory.cpp
      incident_net_array.cpp
      contraction_tree.cpp
      dynamic_graph.cpp
      dynamic_graph_factory.cpp
      dynamic_adjacency_array.cpp
      fixed_vertex_support.cpp
      graph.cpp)

target_sources(MtKaHyPar-ToolsSources INTERFACE ${ToolsDatastructureSources})
