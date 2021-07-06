target_sources(mt_kahypar_fast_tests PRIVATE
        static_hypergraph_test.cc
        dynamic_partitioned_hypergraph_test.cc
        partitioned_hypergraph_test.cc
        delta_partitioned_hypergraph_test.cc
        partitioned_hypergraph_smoke_test.cc
        gain_updates_test.cc
        graph_test.cc
        connectivity_set_test.cc
        priority_queue_test.cc
        sparsifier_hypergraph_test.cc
        array_test.cc
        sparse_map_test.cc
        pin_count_in_part_test.cc
)

set (StrongTestsDatastructuresSources
        dynamic_hypergraph_test.cc
        dynamic_partitioned_hypergraph_test.cc
        partitioned_hypergraph_test.cc
        delta_partitioned_hypergraph_test.cc
        partitioned_hypergraph_smoke_test.cc
        contraction_tree_test.cc
        gain_updates_test.cc
        incident_net_array_test.cc
        nlevel_smoke_test.cc
        smoke_test_common.cc
        )

set(StrongAsyncTestsDatastructuresSources
        uncontraction_group_tree_test.cpp
        group_pool_test.cc
        async_nlevel_smoke_test.cc
        lock_manager_test.cc
        async_common_test.cc
        dynamic_hypergraph_test.cc
        partitioned_hypergraph_test.cc
        delta_partitioned_hypergraph_test.cc
        partitioned_hypergraph_smoke_test.cc
        contraction_tree_test.cc
        gain_updates_test.cc
        incident_net_array_test.cc
        smoke_test_common.cc
        iterable_bit_set_test.cc
        depth_priority_queue_test.cc
        )

target_sources(mt_kahypar_strong_tests PRIVATE ${StrongTestsDatastructuresSources})
target_sources(mt_kahypar_strong_async_tests PRIVATE ${StrongAsyncTestsDatastructuresSources})

