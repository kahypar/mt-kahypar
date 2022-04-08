include(GNUInstallDirs)

# Library MT-KaHyPar-D
add_library(mtkahypard SHARED libmtkahypar.cc)
target_link_libraries(mtkahypard ${Boost_LIBRARIES})
target_compile_definitions(mtkahypard PUBLIC MT_KAHYPAR_LIBRARY_MODE)

# Library MT-KaHyPar-Q
add_library(mtkahyparq SHARED libmtkahypar.cc)
target_link_libraries(mtkahyparq ${Boost_LIBRARIES})
target_compile_definitions(mtkahyparq PUBLIC MT_KAHYPAR_LIBRARY_MODE)
target_compile_definitions(mtkahyparq PUBLIC USE_STRONG_PARTITIONER)

# Library MT-KaHyPar-Graph
add_library(mtkahypargraph SHARED libmtkahypar.cc)
target_link_libraries(mtkahypargraph ${Boost_LIBRARIES})
target_compile_definitions(mtkahypargraph PUBLIC MT_KAHYPAR_LIBRARY_MODE)
target_compile_definitions(mtkahypargraph PUBLIC USE_GRAPH_PARTITIONER)

# Library MT-KaHyPar-GQ
add_library(mtkahypargq SHARED libmtkahypar.cc)
target_link_libraries(mtkahypargq ${Boost_LIBRARIES})
target_compile_definitions(mtkahypargq PUBLIC MT_KAHYPAR_LIBRARY_MODE)
target_compile_definitions(mtkahypargq PUBLIC USE_GRAPH_PARTITIONER USE_STRONG_PARTITIONER)

set(TARGETS_WANTING_ALL_SOURCES ${TARGETS_WANTING_ALL_SOURCES} mtkahypard mtkahyparq mtkahypargraph mtkahypargq PARENT_SCOPE)

set_target_properties(mtkahypard PROPERTIES
    PUBLIC_HEADER ../include/libmtkahypar.h)
set_target_properties(mtkahyparq PROPERTIES
    PUBLIC_HEADER ../include/libmtkahypar.h)
set_target_properties(mtkahypargraph PROPERTIES
    PUBLIC_HEADER ../include/libmtkahypar.h)
set_target_properties(mtkahypargq PROPERTIES
    PUBLIC_HEADER ../include/libmtkahypar.h)

target_include_directories(mtkahypard PRIVATE ../include)
target_include_directories(mtkahyparq PRIVATE ../include)
target_include_directories(mtkahypargraph PRIVATE ../include)
target_include_directories(mtkahypargq PRIVATE ../include)

configure_file(libmtkahypard.pc.in libmtkahypard.pc @ONLY)
configure_file(libmtkahyparq.pc.in libmtkahyparq.pc @ONLY)
configure_file(libmtkahypargraph.pc.in libmtkahypargraph.pc @ONLY)
configure_file(libmtkahypargq.pc.in libmtkahypargq.pc @ONLY)

install(TARGETS mtkahypard mtkahyparq mtkahypargraph mtkahypargq
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(FILES ${CMAKE_BINARY_DIR}/lib/libmtkahypard.pc
              ${CMAKE_BINARY_DIR}/lib/libmtkahyparq.pc
              ${CMAKE_BINARY_DIR}/lib/libmtkahypargraph.pc
              ${CMAKE_BINARY_DIR}/lib/libmtkahypargq.pc
        DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)

configure_file(cmake_uninstall.cmake.in cmake_uninstall.cmake IMMEDIATE @ONLY)

add_custom_target(uninstall-mtkahypar "${CMAKE_COMMAND}" -P cmake_uninstall.cmake)

add_custom_target(install.mtkahypar
    ${CMAKE_COMMAND}
    -DBUILD_TYPE=${CMAKE_BUILD_TYPE}
    -P ${CMAKE_BINARY_DIR}/cmake_install.cmake)
ADD_DEPENDENCIES(install.mtkahypar mtkahypard mtkahyparq mtkahypargraph mtkahypargq)