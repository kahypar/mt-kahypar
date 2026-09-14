if(KAHYPAR_ENABLE_DETERMINISTIC_FEATURES)
  target_sources(mtkahypar_tests PRIVATE
        determinism_test.cc
        )
endif()
