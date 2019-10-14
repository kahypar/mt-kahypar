add_subdirectory(preprocessing)
add_subdirectory(coarsening)
add_subdirectory(refinement)

file(COPY test_instances DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
add_gmock_test(partitioner_test partitioner_test.cc)