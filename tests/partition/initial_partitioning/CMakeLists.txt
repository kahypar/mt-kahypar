target_sources(mt_kahypar_fast_tests PRIVATE
        flat_initial_partitioner_test.cc
        initial_partitioner_test.cc
        initial_partitioning_data_container_test.cc
        pool_initial_partitioner_test.cc
        )

set(StrongTestsIPSources
        flat_initial_partitioner_test.cc
        initial_partitioner_test.cc
        initial_partitioning_data_container_test.cc
        pool_initial_partitioner_test.cc
        )

target_sources(mt_kahypar_strong_tests PRIVATE ${StrongTestsIPSources})
target_sources(mt_kahypar_strong_async_tests PRIVATE ${StrongTestsIPSources})