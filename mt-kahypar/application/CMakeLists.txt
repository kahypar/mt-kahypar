configure_file(mt_kahypar.cc.in ${PROJECT_BINARY_DIR}/mt-kahypar/application/mt_kahypar.cc)
add_executable(MtKaHyParWrapper ${PROJECT_BINARY_DIR}/mt-kahypar/application/mt_kahypar.cc)
set_target_properties(MtKaHyParWrapper PROPERTIES OUTPUT_NAME "MtKaHyPar")
target_link_libraries(MtKaHyParWrapper ${Boost_LIBRARIES})
target_link_libraries(MtKaHyParWrapper pthread)
set_property(TARGET MtKaHyParWrapper PROPERTY CXX_STANDARD 17)
set_property(TARGET MtKaHyParWrapper PROPERTY CXX_STANDARD_REQUIRED ON)

add_executable(MtKaHyParDefault kahypar.cc)
target_link_libraries(MtKaHyParDefault ${Boost_LIBRARIES})
target_link_libraries(MtKaHyParDefault pthread)
set_property(TARGET MtKaHyParDefault PROPERTY CXX_STANDARD 17)
set_property(TARGET MtKaHyParDefault PROPERTY CXX_STANDARD_REQUIRED ON)
add_custom_command(TARGET MtKaHyParDefault PRE_BUILD
    COMMAND ${CMAKE_COMMAND}
      -E echo_append "${MT_KAHYPAR_VERSION_GIT_REFSPEC} at sha ${MT_KAHYPAR_VERSION_GIT_SHA1}" >
      ${PROJECT_BINARY_DIR}/mt-kahypar/application/git_mt_kahypar_d.txt)

add_executable(MtKaHyParQuality kahypar.cc)
target_link_libraries(MtKaHyParQuality ${Boost_LIBRARIES})
target_link_libraries(MtKaHyParQuality pthread)
set_property(TARGET MtKaHyParQuality PROPERTY CXX_STANDARD 17)
set_property(TARGET MtKaHyParQuality PROPERTY CXX_STANDARD_REQUIRED ON)
target_compile_definitions(MtKaHyParQuality PUBLIC USE_STRONG_PARTITIONER)
add_custom_command(TARGET MtKaHyParQuality PRE_BUILD
    COMMAND ${CMAKE_COMMAND}
      -E echo_append "${MT_KAHYPAR_VERSION_GIT_REFSPEC} at sha ${MT_KAHYPAR_VERSION_GIT_SHA1}" >
      ${PROJECT_BINARY_DIR}/mt-kahypar/application/git_mt_kahypar_q.txt)

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

add_executable(MtKaHyParGraphQuality kahypar.cc)
target_link_libraries(MtKaHyParGraphQuality ${Boost_LIBRARIES})
target_link_libraries(MtKaHyParGraphQuality pthread)
set_property(TARGET MtKaHyParGraphQuality PROPERTY CXX_STANDARD 17)
set_property(TARGET MtKaHyParGraphQuality PROPERTY CXX_STANDARD_REQUIRED ON)
target_compile_definitions(MtKaHyParGraphQuality PUBLIC USE_GRAPH_PARTITIONER USE_STRONG_PARTITIONER)
add_custom_command(TARGET MtKaHyParGraphQuality PRE_BUILD
    COMMAND ${CMAKE_COMMAND}
      -E echo_append "${MT_KAHYPAR_VERSION_GIT_REFSPEC} at sha ${MT_KAHYPAR_VERSION_GIT_SHA1}" >
      ${PROJECT_BINARY_DIR}/mt-kahypar/application/git_mt_kahypar_gq.txt)

if(ENABLE_PROFILE MATCHES ON)
  target_link_libraries(MtKaHyParWrapper ${PROFILE_FLAGS})
  target_link_libraries(MtKaHyParDefault ${PROFILE_FLAGS})
  target_link_libraries(MtKaHyParQuality ${PROFILE_FLAGS})
  target_link_libraries(MtKaHyParGraph ${PROFILE_FLAGS})
  target_link_libraries(MtKaHyParGraphQuality ${PROFILE_FLAGS})
endif()


set(TARGETS_WANTING_ALL_SOURCES ${TARGETS_WANTING_ALL_SOURCES} MtKaHyParWrapper MtKaHyParDefault MtKaHyParQuality MtKaHyParGraph MtKaHyParGraphQuality PARENT_SCOPE)