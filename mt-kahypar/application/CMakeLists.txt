configure_file(mt_kahypar.cc.in ${PROJECT_BINARY_DIR}/mt-kahypar/application/mt_kahypar.cc)
add_executable(MtKaHyPar ${PROJECT_BINARY_DIR}/mt-kahypar/application/mt_kahypar.cc)
target_link_libraries(MtKaHyPar ${Boost_LIBRARIES})
target_link_libraries(MtKaHyPar pthread)
set_property(TARGET MtKaHyPar PROPERTY CXX_STANDARD 17)
set_property(TARGET MtKaHyPar PROPERTY CXX_STANDARD_REQUIRED ON)

if(ENABLE_PROFILE MATCHES ON)
  target_link_libraries(MtKaHyParDefault ${PROFILE_FLAGS})
endif()

set(TARGETS_WANTING_ALL_SOURCES ${TARGETS_WANTING_ALL_SOURCES} MtKaHyPar PARENT_SCOPE)