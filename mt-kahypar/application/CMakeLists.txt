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

add_executable(MtKaHyParFastGraph kahypar.cc)
target_link_libraries(MtKaHyParFastGraph ${Boost_LIBRARIES})
target_link_libraries(MtKaHyParFastGraph pthread)
set_property(TARGET MtKaHyParFastGraph PROPERTY CXX_STANDARD 17)
set_property(TARGET MtKaHyParFastGraph PROPERTY CXX_STANDARD_REQUIRED ON)
target_compile_definitions(MtKaHyParFastGraph PUBLIC USE_GRAPH_STRUCTURE)

if(ENABLE_PROFILE MATCHES ON)
  target_link_libraries(MtKaHyParFast ${PROFILE_FLAGS})
  target_link_libraries(MtKaHyParStrong ${PROFILE_FLAGS})
  target_link_libraries(MtKaHyParFastGraph ${PROFILE_FLAGS})
endif()


set(TARGETS_WANTING_ALL_SOURCES ${TARGETS_WANTING_ALL_SOURCES} MtKaHyParFast MtKaHyParStrong MtKaHyParFastGraph PARENT_SCOPE)