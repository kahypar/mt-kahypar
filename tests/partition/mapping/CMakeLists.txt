if ( KAHYPAR_ENABLE_STEINER_TREE_METRIC )
  target_sources(mtkahypar_tests PRIVATE
          target_graph_test.cc
          set_enumerator_test.cc
          initial_mapping_test.cc
          )
endif()
