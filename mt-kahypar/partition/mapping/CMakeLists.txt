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

set(ToolsMappingSources
        target_graph.cpp
        all_pair_shortest_path.cpp
        steiner_tree.cpp
        )

foreach(modtarget IN LISTS TOOLS_TARGETS)
    target_sources(${modtarget} PRIVATE ${ToolsMappingSources})
endforeach()

set(OneToOneMappingSources
      greedy_mapping.cpp
      initial_mapping.cpp
      kerninghan_lin.cpp
      )
target_sources(OneToOneMapping PRIVATE ${OneToOneMappingSources})
