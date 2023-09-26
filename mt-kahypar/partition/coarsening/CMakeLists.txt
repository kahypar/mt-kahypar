set(MultilevelCoarseningSources
        multilevel/deterministic_multilevel_coarsener.cpp
        multilevel/multilevel_uncoarsener.cpp)

set(NLevelCoarseningSources
        nlevel/nlevel_uncoarsener.cpp)

foreach(modtarget IN LISTS PARTITIONING_SUITE_TARGETS)
    target_sources(${modtarget} PRIVATE ${MultilevelCoarseningSources})
    target_sources(${modtarget} PRIVATE ${NLevelCoarseningSources})
endforeach()