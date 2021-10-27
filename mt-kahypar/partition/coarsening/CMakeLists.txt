set(CoarseningSources
        deterministic_multilevel_coarsener.cpp
        multilevel_uncoarsener.cpp
        nlevel_uncoarsener.cpp
        )

foreach(modtarget IN LISTS TARGETS_WANTING_ALL_SOURCES)
    target_sources(${modtarget} PRIVATE ${CoarseningSources})
endforeach()
