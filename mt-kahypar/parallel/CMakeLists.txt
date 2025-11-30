set(ParallelSources
        thread_management.cpp)

target_sources(MtKaHyPar-Sources INTERFACE ${ParallelSources})
target_sources(MtKaHyPar-ToolsSources INTERFACE ${ParallelSources})
