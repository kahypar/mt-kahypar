set(PreprocessingSources
        community_detection/parallel_louvain.cpp
        community_detection/local_moving_modularity.cpp)

target_sources(KaHyPar PRIVATE ${PreprocessingSources})
target_sources(mt_kahypar_tests PRIVATE ${PreprocessingSources})