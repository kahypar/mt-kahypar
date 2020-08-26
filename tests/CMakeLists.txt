add_executable(mt_kahypar_tests run_tests.cpp)
target_link_libraries(mt_kahypar_tests gtest gtest_main ${CMAKE_THREAD_LIBS_INIT})
set_property(TARGET mt_kahypar_tests PROPERTY CXX_STANDARD 17)
set_property(TARGET mt_kahypar_tests PROPERTY CXX_STANDARD_REQUIRED ON)

add_subdirectory(datastructures)
add_subdirectory(end_to_end)
add_subdirectory(io)
add_subdirectory(parallel)
add_subdirectory(partition)

set(TARGETS_WANTING_ALL_SOURCES ${TARGETS_WANTING_ALL_SOURCES} mt_kahypar_tests PARENT_SCOPE)
