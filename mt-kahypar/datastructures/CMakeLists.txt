set(DatastructureSources static_hypergraph_factory.cpp static_hypergraph.cpp graph.cpp)

target_sources(KaHyPar PRIVATE ${DatastructureSources})
target_sources(mt_kahypar_tests PRIVATE ${DatastructureSources})       # TODO make list collecting targets that want everything?