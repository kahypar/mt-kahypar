configure_file(mt_kahypar.cc.in ${PROJECT_BINARY_DIR}/mt-kahypar/application/mt_kahypar.cc)
add_executable(MtKaHyPar ${PROJECT_BINARY_DIR}/mt-kahypar/application/mt_kahypar.cc)
target_link_libraries(MtKaHyPar ${Boost_LIBRARIES})
target_link_libraries(MtKaHyPar pthread)
set_property(TARGET MtKaHyPar PROPERTY CXX_STANDARD 17)
set_property(TARGET MtKaHyPar PROPERTY CXX_STANDARD_REQUIRED ON)

add_executable(MtKaHyParFast kahypar.cc)
target_link_libraries(MtKaHyParFast ${Boost_LIBRARIES})
target_link_libraries(MtKaHyParFast pthread)
set_property(TARGET MtKaHyParFast PROPERTY CXX_STANDARD 17)
set_property(TARGET MtKaHyParFast PROPERTY CXX_STANDARD_REQUIRED ON)
add_custom_command(TARGET MtKaHyParFast PRE_BUILD
    COMMAND ${CMAKE_COMMAND}
      -E echo_append "${MT_KAHYPAR_VERSION_GIT_REFSPEC} at sha ${MT_KAHYPAR_VERSION_GIT_SHA1}" >
      ${PROJECT_BINARY_DIR}/mt-kahypar/application/git_mt_kahypar_fast.txt)

add_executable(MtKaHyParStrong kahypar.cc)
target_link_libraries(MtKaHyParStrong ${Boost_LIBRARIES})
target_link_libraries(MtKaHyParStrong pthread)
set_property(TARGET MtKaHyParStrong PROPERTY CXX_STANDARD 17)
set_property(TARGET MtKaHyParStrong PROPERTY CXX_STANDARD_REQUIRED ON)
target_compile_definitions(MtKaHyParStrong PUBLIC USE_STRONG_PARTITIONER)
add_custom_command(TARGET MtKaHyParStrong PRE_BUILD
    COMMAND ${CMAKE_COMMAND}
      -E echo_append "${MT_KAHYPAR_VERSION_GIT_REFSPEC} at sha ${MT_KAHYPAR_VERSION_GIT_SHA1}" >
      ${PROJECT_BINARY_DIR}/mt-kahypar/application/git_mt_kahypar_strong.txt)

add_executable(MtKaHyParGraph kahypar.cc)
target_link_libraries(MtKaHyParGraph ${Boost_LIBRARIES})
target_link_libraries(MtKaHyParGraph pthread)
set_property(TARGET MtKaHyParGraph PROPERTY CXX_STANDARD 17)
set_property(TARGET MtKaHyParGraph PROPERTY CXX_STANDARD_REQUIRED ON)
target_compile_definitions(MtKaHyParGraph PUBLIC USE_GRAPH_PARTITIONER)
add_custom_command(TARGET MtKaHyParGraph PRE_BUILD
    COMMAND ${CMAKE_COMMAND}
      -E echo_append "${MT_KAHYPAR_VERSION_GIT_REFSPEC} at sha ${MT_KAHYPAR_VERSION_GIT_SHA1}" >
      ${PROJECT_BINARY_DIR}/mt-kahypar/application/git_mt_kahypar_graph.txt)

if(ENABLE_PROFILE MATCHES ON)
  target_link_libraries(MtKaHyParFast ${PROFILE_FLAGS})
  target_link_libraries(MtKaHyParStrong ${PROFILE_FLAGS})
  target_link_libraries(MtKaHyParGraph ${PROFILE_FLAGS})
endif()


set(TARGETS_WANTING_ALL_SOURCES ${TARGETS_WANTING_ALL_SOURCES} MtKaHyPar MtKaHyParFast MtKaHyParStrong MtKaHyParGraph PARENT_SCOPE)