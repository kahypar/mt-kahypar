add_gmock_test(graph_test graph_test.cc)
add_gmock_test(static_hypergraph_test static_hypergraph_test.cc)
add_gmock_test(partitioned_hypergraph_test partitioned_hypergraph_test.cc)
add_gmock_test(delta_partitioned_hypergraph_test delta_partitioned_hypergraph_test.cc)
add_gmock_test(partitioned_hypergraph_smoke_test partitioned_hypergraph_smoke_test.cc)
add_gmock_test(connectivity_set_test connectivity_set_test.cc)
add_gmock_test(priority_queue_test priority_queue_test.cc)
add_gmock_test(gain_updates_test gain_updates_test.cc)
add_gmock_test(sparsifier_hypergraph_test sparsifier_hypergraph_test.cc)
add_gmock_test(array_test array_test.cc)
add_gmock_test(pin_count_in_part_test pin_count_in_part_test.cc)
add_gmock_test(sparse_map_test sparse_map_test.cc)

target_sources(mt_kahypar_tests PRIVATE
        static_hypergraph_test.cc
        connectivity_set_test.cc
        array_test.cc
        priority_queue_test.cc
        )
