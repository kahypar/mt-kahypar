target_sources(mt_kahypar_tests PRIVATE
        flat_initial_partitioner_test.cc
        initial_partitioning_data_container_test.cc
        pool_initial_partitioner_test.cc
        )

if(NOT MT_KAHYPAR_DISABLE_BOOST)
        target_sources(mt_kahypar_tests PRIVATE initial_partitioner_test.cc)
endif()
