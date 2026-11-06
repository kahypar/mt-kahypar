add_gmock_test(mtkahypar_interface_test interface_test.cc)
target_link_libraries(mtkahypar_interface_test MtKaHyPar-LibraryBuildSources MtKaHyPar-Test mtkahypar)
# The test reads its input files from the source tree, such that it can be executed from any directory.
target_compile_definitions(mtkahypar_interface_test PRIVATE
  "MT_KAHYPAR_INTERFACE_TEST_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}\"")
