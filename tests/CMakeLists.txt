# Contains tests specific to our fast partitioner
add_executable(mt_kahypar_multilevel_tests run_tests.cpp)
target_link_libraries(mt_kahypar_multilevel_tests gtest gtest_main ${CMAKE_THREAD_LIBS_INIT})
set_property(TARGET mt_kahypar_multilevel_tests PROPERTY CXX_STANDARD 17)
set_property(TARGET mt_kahypar_multilevel_tests PROPERTY CXX_STANDARD_REQUIRED ON)
target_link_libraries(mt_kahypar_multilevel_tests ${Boost_LIBRARIES})

# Contains tests specific to our strong partitioner
add_executable(mt_kahypar_nlevel_tests run_tests.cpp)
target_link_libraries(mt_kahypar_nlevel_tests gtest gtest_main ${CMAKE_THREAD_LIBS_INIT})
set_property(TARGET mt_kahypar_nlevel_tests PROPERTY CXX_STANDARD 17)
set_property(TARGET mt_kahypar_nlevel_tests PROPERTY CXX_STANDARD_REQUIRED ON)
target_compile_definitions(mt_kahypar_nlevel_tests PUBLIC USE_STRONG_PARTITIONER)
target_link_libraries(mt_kahypar_nlevel_tests ${Boost_LIBRARIES})

# Contains tests specific to our fast partitioner for graphs
add_executable(mt_kahypar_graph_tests run_tests.cpp)
target_link_libraries(mt_kahypar_graph_tests gtest gtest_main ${CMAKE_THREAD_LIBS_INIT})
set_property(TARGET mt_kahypar_graph_tests PROPERTY CXX_STANDARD 17)
set_property(TARGET mt_kahypar_graph_tests PROPERTY CXX_STANDARD_REQUIRED ON)
target_compile_definitions(mt_kahypar_graph_tests PUBLIC USE_GRAPH_PARTITIONER)
target_link_libraries(mt_kahypar_graph_tests ${Boost_LIBRARIES})

add_subdirectory(datastructures)
add_subdirectory(interface)
add_subdirectory(io)
add_subdirectory(parallel)
add_subdirectory(partition)

set(TARGETS_WANTING_ALL_SOURCES ${TARGETS_WANTING_ALL_SOURCES} mt_kahypar_multilevel_tests mt_kahypar_nlevel_tests mt_kahypar_graph_tests PARENT_SCOPE)
