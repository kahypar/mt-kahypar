# Contains tests specific to our fast partitioner
add_executable(mt_kahypar_fast_tests run_tests.cpp)
target_link_libraries(mt_kahypar_fast_tests gtest gtest_main ${CMAKE_THREAD_LIBS_INIT})
set_property(TARGET mt_kahypar_fast_tests PROPERTY CXX_STANDARD 17)
set_property(TARGET mt_kahypar_fast_tests PROPERTY CXX_STANDARD_REQUIRED ON)
target_link_libraries(mt_kahypar_fast_tests ${Boost_LIBRARIES})

# Contains tests specific to our strong partitioner
add_executable(mt_kahypar_strong_tests run_tests.cpp)
target_link_libraries(mt_kahypar_strong_tests gtest gtest_main ${CMAKE_THREAD_LIBS_INIT})
set_property(TARGET mt_kahypar_strong_tests PROPERTY CXX_STANDARD 17)
set_property(TARGET mt_kahypar_strong_tests PROPERTY CXX_STANDARD_REQUIRED ON)
target_compile_definitions(mt_kahypar_strong_tests PUBLIC USE_STRONG_PARTITIONER)
target_link_libraries(mt_kahypar_strong_tests ${Boost_LIBRARIES})

add_subdirectory(datastructures)
add_subdirectory(end_to_end)
add_subdirectory(io)
add_subdirectory(parallel)
add_subdirectory(partition)

set(TARGETS_WANTING_ALL_SOURCES ${TARGETS_WANTING_ALL_SOURCES} mt_kahypar_fast_tests mt_kahypar_strong_tests PARENT_SCOPE)
