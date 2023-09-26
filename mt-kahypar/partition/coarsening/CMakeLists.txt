set(MultilevelCoarseningSources
        multilevel/deterministic_multilevel_coarsener.cpp
        multilevel/multilevel_uncoarsener.cpp
        multilevel/concurrent_clustering_data.cpp)

set(NLevelCoarseningSources
        nlevel/nlevel_uncoarsener.cpp)

target_sources(MtKaHyPar-Sources INTERFACE ${MultilevelCoarseningSources})
target_sources(MtKaHyPar-Sources INTERFACE ${NLevelCoarseningSources})
