target_sources(mtkahypar_tests PRIVATE
        flat_initial_partitioner_test.cc
        initial_partitioning_data_container_test.cc
        )

if(NOT KAHYPAR_DISABLE_HWLOC)
        target_sources(mtkahypar_tests PRIVATE
                initial_partitioner_test.cc
                pool_initial_partitioner_test.cc)
endif()
