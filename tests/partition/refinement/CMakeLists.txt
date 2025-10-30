target_sources(mtkahypar_tests PRIVATE
         bipartitioning_gain_policy_test.cc
         deterministic_jet_refiner_test.cc
         gain_policy_test.cc
         rollback_test.cc
         advanced_rebalancer_test.cc
         deterministic_rebalancer_test.cc
         twoway_fm_refiner_test.cc
         gain_cache_test.cc
         multitry_fm_test.cc
         fm_strategy_test.cc
         )

if(NOT KAHYPAR_DISABLE_HWLOC)
    target_sources(mtkahypar_tests PRIVATE label_propagation_refiner_test.cc)
endif()
if(NOT KAHYPAR_MINIMAL_COMPILATION)
    target_sources(mtkahypar_tests PRIVATE
            problem_construction_test.cc
            scheduler_test.cc
            rebalance_test.cc
            flow_construction_test.cc
            )
endif()
