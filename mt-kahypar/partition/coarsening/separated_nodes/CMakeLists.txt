set(SeparatedNodesSources
        snodes_coarsening_pass.cpp
        )

foreach(modtarget IN LISTS TARGETS_WANTING_ALL_SOURCES)
    target_sources(${modtarget} PRIVATE ${SeparatedNodesSources})
endforeach()
