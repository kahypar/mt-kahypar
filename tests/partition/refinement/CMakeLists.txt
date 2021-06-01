target_sources(mt_kahypar_fast_tests PRIVATE
        gain_policy_test.cc
        label_propagation_refiner_test.cc
        rollback_test.cc
        rebalance_test.cc
        twoway_fm_refiner_test.cc
        gain_test.cc
        multitry_fm_test.cc
        fm_strategy_test.cc
        )

set(StrongTestsRefinementSources
        gain_policy_test.cc
        label_propagation_refiner_test.cc
        rollback_test.cc
        rebalance_test.cc
        twoway_fm_refiner_test.cc
        gain_test.cc
        multitry_fm_test.cc
        fm_strategy_test.cc
        )

set(StrongAsyncTestsRefinementSources
        async_lp_refiner_test.cc
        )

target_sources(mt_kahypar_strong_tests PRIVATE ${StrongTestsRefinementSources})
target_sources(mt_kahypar_strong_async_tests PRIVATE ${StrongTestsRefinementSources} ${StrongAsyncTestsRefinementSources})
