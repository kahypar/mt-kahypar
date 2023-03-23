set(UtilSources
      memory_tree.cpp
    )

foreach(modtarget IN LISTS TARGETS_WANTING_ALL_SOURCES)
    target_sources(${modtarget} PRIVATE ${UtilSources})
endforeach()