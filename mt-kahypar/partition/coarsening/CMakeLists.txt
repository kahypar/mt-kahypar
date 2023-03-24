set(CoarseningSources
        deterministic_multilevel_coarsener.cpp
        multilevel_uncoarsener.cpp
        nlevel_uncoarsener.cpp
        )

foreach(modtarget IN LISTS PARTITIONING_SUITE_TARGETS)
    target_sources(${modtarget} PRIVATE ${CoarseningSources})
endforeach()
