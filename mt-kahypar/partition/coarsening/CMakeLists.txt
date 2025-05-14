set(MultilevelCoarseningSources
        deterministic_multilevel_coarsener.cpp
        multilevel_uncoarsener.cpp
        experimental_coarsener.cpp)

set(NLevelCoarseningSources
        nlevel_uncoarsener.cpp)

target_sources(MtKaHyPar-Sources INTERFACE ${MultilevelCoarseningSources})
target_sources(MtKaHyPar-Sources INTERFACE ${NLevelCoarseningSources})
