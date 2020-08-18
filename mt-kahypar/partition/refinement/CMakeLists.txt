set(RefinementSources
        fm/multitry_kway_fm.cpp
        fm/localized_kway_fm_core.cpp
        fm/global_rollback.cpp
        label_propagation/label_propagation_refiner.cpp
        rebalancing/rebalancer.cpp
        )

target_sources(KaHyPar PRIVATE ${RefinementSources})
target_sources(mt_kahypar_tests PRIVATE ${RefinementSources})
