target_sources(mt_kahypar_multilevel_tests PRIVATE
        refinement_adapter_test.cc
        problem_construction_test.cc
        scheduler_test.cc
        quotient_graph_test.cc
        gain_policy_test.cc
        label_propagation_refiner_test.cc
        rollback_test.cc
        rebalance_test.cc
        twoway_fm_refiner_test.cc
        gain_test.cc
        multitry_fm_test.cc
        fm_strategy_test.cc
        flow_construction_test.cc
        )

target_sources(mt_kahypar_nlevel_tests PRIVATE
        refinement_adapter_test.cc
        problem_construction_test.cc
        scheduler_test.cc
        quotient_graph_test.cc
        gain_policy_test.cc
        label_propagation_refiner_test.cc
        rollback_test.cc
        rebalance_test.cc
        twoway_fm_refiner_test.cc
        gain_test.cc
        multitry_fm_test.cc
        fm_strategy_test.cc
        flow_construction_test.cc
        )
