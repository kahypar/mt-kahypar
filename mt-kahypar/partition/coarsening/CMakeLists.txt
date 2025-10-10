set(MultilevelCoarseningSources
        multilevel/multilevel_uncoarsener.cpp
        multilevel/concurrent_clustering_data.cpp)

set(DeterministicCoarseningSources
        multilevel/deterministic_multilevel_coarsener.cpp)

set(NLevelCoarseningSources
        nlevel/nlevel_uncoarsener.cpp)

target_sources(MtKaHyPar-Sources INTERFACE ${MultilevelCoarseningSources})

if (NOT KAHYPAR_MINIMAL_COMPILATION)
        target_sources(MtKaHyPar-Sources INTERFACE ${DeterministicCoarseningSources})
        target_sources(MtKaHyPar-Sources INTERFACE ${NLevelCoarseningSources})
endif()
