target_sources(mt_kahypar_tests PRIVATE
        flat_initial_partitioner_test.cc
        initial_partitioner_test.cc
        initial_partitioning_data_container_test.cc
        pool_initial_partitioner_test.cc
        )

target_link_libraries(mt_kahypar_tests ${Boost_LIBRARIES})

add_gmock_test(flat_initial_partitioner_test flat_initial_partitioner_test.cc)
add_gmock_test(initial_partitioner_test initial_partitioner_test.cc)
add_gmock_test(initial_partitioning_data_container_test initial_partitioning_data_container_test.cc)
target_link_libraries(initial_partitioner_test ${Boost_LIBRARIES})
add_gmock_test(pool_initial_partitioner_test pool_initial_partitioner_test.cc)