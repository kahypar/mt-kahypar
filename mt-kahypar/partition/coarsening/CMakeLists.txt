set(CoarseningSources extended_clustering.cpp uncoarsening.cpp)

target_sources(KaHyPar PRIVATE ${CoarseningSources})
target_sources(mt_kahypar_tests PRIVATE ${CoarseningSources})