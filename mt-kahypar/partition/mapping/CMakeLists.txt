set(MappingSources
        target_graph.cpp
        all_pair_shortest_path.cpp
        steiner_tree.cpp
        greedy_mapping.cpp
        initial_mapping.cpp
        kerninghan_lin.cpp
        )

foreach(modtarget IN LISTS PARTITIONING_SUITE_TARGETS)
    target_sources(${modtarget} PRIVATE ${MappingSources})
endforeach()

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