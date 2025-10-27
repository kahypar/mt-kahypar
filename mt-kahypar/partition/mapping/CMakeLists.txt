set(MappingSources
        target_graph.cpp
        all_pair_shortest_path.cpp
        steiner_tree.cpp
        greedy_mapping.cpp
        initial_mapping.cpp
        kerninghan_lin.cpp
        )

if ( KAHYPAR_ENABLE_STEINER_TREE_METRIC )
  target_sources(MtKaHyPar-Sources INTERFACE ${MappingSources})
else ()
  target_sources(MtKaHyPar-LibraryBuildSources PRIVATE ${MappingSources})
endif()

if(NOT KAHYPAR_MINIMAL_COMPILATION)
  target_sources(MtKaHyPar-ToolsSources INTERFACE ${MappingSources})
endif()
