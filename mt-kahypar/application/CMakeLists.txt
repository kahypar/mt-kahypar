add_executable(MtKaHyParFast kahypar.cc)
target_link_libraries(MtKaHyParFast ${Boost_LIBRARIES})
target_link_libraries(MtKaHyParFast pthread)
set_property(TARGET MtKaHyParFast PROPERTY CXX_STANDARD 17)
set_property(TARGET MtKaHyParFast PROPERTY CXX_STANDARD_REQUIRED ON)

add_executable(MtKaHyParStrong kahypar.cc)
target_link_libraries(MtKaHyParStrong ${Boost_LIBRARIES})
target_link_libraries(MtKaHyParStrong pthread)
set_property(TARGET MtKaHyParStrong PROPERTY CXX_STANDARD 17)
set_property(TARGET MtKaHyParStrong PROPERTY CXX_STANDARD_REQUIRED ON)
target_compile_definitions(MtKaHyParStrong PUBLIC USE_STRONG_PARTITIONER)

add_executable(MtKaHyParGraph kahypar.cc)
target_link_libraries(MtKaHyParGraph ${Boost_LIBRARIES})
target_link_libraries(MtKaHyParGraph pthread)
set_property(TARGET MtKaHyParGraph PROPERTY CXX_STANDARD 17)
set_property(TARGET MtKaHyParGraph PROPERTY CXX_STANDARD_REQUIRED ON)
target_compile_definitions(MtKaHyParGraph PUBLIC USE_GRAPH_STRUCTURE)

if(ENABLE_PROFILE MATCHES ON)
  target_link_libraries(MtKaHyParFast ${PROFILE_FLAGS})
  target_link_libraries(MtKaHyParStrong ${PROFILE_FLAGS})
  target_link_libraries(MtKaHyParGraph ${PROFILE_FLAGS})
endif()


set(TARGETS_WANTING_ALL_SOURCES ${TARGETS_WANTING_ALL_SOURCES} MtKaHyParFast MtKaHyParStrong MtKaHyParGraph PARENT_SCOPE)