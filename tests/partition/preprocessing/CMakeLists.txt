target_sources(mt_kahypar_tests PRIVATE
        louvain_test.cc
        )

add_gmock_test(similiar_net_combiner_test similiar_net_combiner_test.cc)