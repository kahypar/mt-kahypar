set(PreprocessingSources
        community_detection/parallel_louvain.cpp
        community_detection/local_moving_modularity.cpp)

target_sources(MtKaHyPar-Sources INTERFACE ${PreprocessingSources})
target_sources(MtKaHyPar-ToolsSources INTERFACE ${PreprocessingSources})