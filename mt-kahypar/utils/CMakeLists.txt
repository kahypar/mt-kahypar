set(UtilSources
      memory_tree.cpp
    )

target_sources(MtKaHyPar-Sources INTERFACE ${UtilSources})
target_sources(MtKaHyPar-ToolsSources INTERFACE ${UtilSources})
