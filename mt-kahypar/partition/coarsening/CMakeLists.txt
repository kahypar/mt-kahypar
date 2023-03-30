set(MultilevelCoarseningSources
        deterministic_multilevel_coarsener.cpp
        multilevel_uncoarsener.cpp)

set(NLevelCoarseningSources
        nlevel_uncoarsener.cpp)

foreach(modtarget IN LISTS PARTITIONING_SUITE_TARGETS)
    target_sources(${modtarget} PRIVATE ${MultilevelCoarseningSources})
endforeach()

if ( KAHYPAR_ENABLE_N_LEVEL_PARTITIONING_FEATURES )
  foreach(modtarget IN LISTS PARTITIONING_SUITE_TARGETS)
    target_sources(${modtarget} PRIVATE ${NLevelCoarseningSources})
  endforeach()
endif()