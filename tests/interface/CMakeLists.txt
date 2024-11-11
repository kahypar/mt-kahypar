file(COPY test_instances DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_gmock_test(interface_test interface_test.cc)
target_link_libraries(interface_test MtKaHyPar-LibraryBuildSources MtKaHyPar-Test mtkahypar)
