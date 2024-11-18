add_executable(MtKaHyPar mt_kahypar_dynamic.cc)
target_link_libraries(MtKaHyPar ${Boost_LIBRARIES})
target_link_libraries(MtKaHyPar TBB::tbb TBB::tbbmalloc_proxy)
target_link_libraries(MtKaHyPar pthread)
set_property(TARGET MtKaHyPar PROPERTY CXX_STANDARD 17)
set_property(TARGET MtKaHyPar PROPERTY CXX_STANDARD_REQUIRED ON)

if(ENABLE_PROFILE MATCHES ON)
  target_link_libraries(MtKaHyPar ${PROFILE_FLAGS})
endif()

set(PARTITIONING_SUITE_TARGETS ${PARTITIONING_SUITE_TARGETS} MtKaHyPar PARENT_SCOPE)