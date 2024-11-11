set(UtilSources
      memory_tree.cpp
    )

target_sources(MtKaHyPar-Sources INTERFACE ${UtilSources})

foreach(modtarget IN LISTS TOOLS_TARGETS)
    target_sources(${modtarget} PRIVATE ${UtilSources})
endforeach()