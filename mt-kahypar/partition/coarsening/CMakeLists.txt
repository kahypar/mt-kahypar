set(MultilevelCoarseningSources
        deterministic_multilevel_coarsener.cpp
        multilevel_uncoarsener.cpp)

set(NLevelCoarseningSources
        nlevel_uncoarsener.cpp)

target_sources(MtKaHyPar-Sources INTERFACE ${MultilevelCoarseningSources})

if(KAHYPAR_ENABLE_HIGHEST_QUALITY_FEATURES)
  target_sources(MtKaHyPar-Sources INTERFACE ${NLevelCoarseningSources})
else()
  target_sources(MtKaHyPar-LibraryBuildSources PRIVATE ${NLevelCoarseningSources})
endif()
